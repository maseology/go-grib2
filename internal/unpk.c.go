package internal

/*
 * unpack grib -- only some formats (code table 5.0) are supported
 *
 * supported: 0 (simple), 4 (ieee), 40 (jpeg), 41(png), 42(aec)
 *
 * input:  sec[]
 *         float data[npnts]
 *
 */
func unpk_grib(sec [][]unsigned_char, data []float) error {

	var packing, bitmap_flag, nbits int
	var ndata unsigned_int
	var mask_pointer []unsigned_char
	var mask unsigned_char
	var ieee, p []unsigned_char
	// float reference, tmp;
	var reference double
	var bin_scale, dec_scale, b double
	uzero := unsigned_int(0)

	packing = code_table_5_0(sec)
	// ndata = (int) GB2_Sec3_npts(sec);
	ndata = GB2_Sec3_npts(sec)
	bitmap_flag = code_table_6_0(sec)

	if bitmap_flag != 0 && bitmap_flag != 254 && bitmap_flag != 255 {
		return fatal_error("unknown bitmap", "")
	}

	if packing == 4 { // ieee
		if sec[5][11] != 1 {
			return fatal_error_i("unpk ieee grib file precision %d not supported", int(sec[5][11]))
		}

		// ieee depacking -- simple no bitmap
		if bitmap_flag == 255 {
			for ii := uzero; ii < ndata; ii++ {
				data[ii] = ieee2flt_nan(sec[7][5+ii*4:])
			}
			return nil
		}
		if bitmap_flag == 0 || bitmap_flag == 254 {
			mask_pointer = sec[6][6:]
			ieee = sec[7][5:]
			ieee_index := 0
			mask = 0
			mask_pointer_index := 0
			for ii := uzero; ii < ndata; ii++ {
				if (ii & 7) == 0 {
					mask = mask_pointer[mask_pointer_index]
					mask_pointer_index++
				}
				if (mask & 128) != 0 {
					data[ii] = ieee2flt_nan(ieee[ieee_index:])
					ieee_index += 4
				} else {
					data[ii] = UNDEFINED
				}
				mask <<= 1
			}
			return nil
		}
		return fatal_error("unknown bitmap", "")
	} else if packing == 0 || packing == 61 { // simple grib1 packing  61 -- log preprocessing

		p = sec[5]
		reference = double(ieee2flt(p[11:]))
		bin_scale = Int_Power(2.0, int2(p[15:]))
		dec_scale = Int_Power(10.0, -int2(p[17:]))
		nbits = int(p[19])
		b = 0.0
		if packing == 61 {
			b = double(ieee2flt(p[20:]))
		}

		if bitmap_flag != 0 && bitmap_flag != 254 && bitmap_flag != 255 {
			return fatal_error("unknown bitmap", "")
		}

		if nbits == 0 {
			tmp := float(reference * dec_scale)
			if packing == 61 {
				tmp = float(exp(double(tmp)) - b)
			} // remove log prescaling
			if bitmap_flag == 255 {
				for ii := uzero; ii < ndata; ii++ {
					data[ii] = tmp
				}
				return nil
			}
			if bitmap_flag == 0 || bitmap_flag == 254 {
				mask_pointer = sec[6][6:]
				mask = 0
				mask_pointer_index := 0
				for ii := uzero; ii < ndata; ii++ {
					if (ii & 7) == 0 {
						mask = mask_pointer[mask_pointer_index]
						mask_pointer_index++
					}
					// data[ii] = (mask & 128) ?  tmp : UNDEFINED;
					if (mask & 128) != 0 {
						data[ii] = tmp
					} else {
						data[ii] = UNDEFINED
					}
					mask <<= 1
				}
				return nil
			}
		}

		// mask_pointer = (bitmap_flag == 255) ? NULL : sec[6] + 6;
		if bitmap_flag == 255 {
			mask_pointer = nil
		} else {
			mask_pointer = sec[6][6:]
		}

		unpk_0(data, sec[7][5:], mask_pointer, nbits, ndata, reference,
			bin_scale, dec_scale)

		if packing == 61 { // remove log prescaling
			// #pragma omp parallel for private(ii) schedule(static)
			for ii := uzero; ii < ndata; ii++ {
				if DEFINED_VAL(data[ii]) {
					data[ii] = float(exp(double(data[ii])) - b)
				}
			}
		}
		return nil
	} else if packing == 2 || packing == 3 { // complex
		return fatal_error("unpk_complex is not supported")
		// TODO: unpk_complex
		// return unpk_complex(sec, data, ndata)
	} else if packing == 40 || packing == 40000 { // jpeg2000

		p = sec[5]
		reference = double(ieee2flt(p[11:]))
		bin_scale = Int_Power(2.0, int2(p[15:]))
		dec_scale = Int_Power(10.0, -int2(p[17:]))
		nbits = int(p[19])
		b = 0.0

		if nbits == 0 {
			tmp := float(reference * dec_scale)
			if bitmap_flag == 255 {
				for ii := uzero; ii < ndata; ii++ {
					data[ii] = tmp
				}
				return nil
			}
			if bitmap_flag == 0 || bitmap_flag == 254 {
				mask_pointer = sec[6][6:]
				ieee = sec[7][5:]
				mask = 0
				mask_pointer_index := 0
				for ii := uzero; ii < ndata; ii++ {
					if (ii & 7) == 0 {
						mask = mask_pointer[mask_pointer_index]
						mask_pointer_index++
					}
					// data[ii] = (mask & 128) ?  tmp : UNDEFINED;
					if (mask & 128) != 0 {
						data[ii] = tmp
					} else {
						data[ii] = UNDEFINED
					}
					mask <<= 1
				}
				return nil
			}
			return fatal_error("unknown bitmap", "")
		}

		// decode jpeg2000
		type img struct{ numcmpts, height, width unsigned_int }
		type opt struct{ buffsize, growable unsigned_int }

		image := img{height: 935, width: 824}
		// opts := opt{buffsize: GB2_Sec7_size(sec) - 5}

		// opts := nil
		// jpcstream=jas_stream_memopen((char *) sec[7][5:], (int) GB2_Sec7_size(sec)-5);
		// image = jpc_decode(jpcstream, opts)
		// if image = nil {
		// 	fatal_error("jpeg2000 decoding", "")
		// }
		// pcmpt = image->cmpts_[0];
		// if (image->numcmpts_ != 1 )
		// return fatal_error("unpk: Found color image.  Grayscale expected","");

		// jas_data=jas_matrix_create(jas_image_height(image), jas_image_width(image));
		// jas_image_readcmpt(image,0,0,0,jas_image_width(image), jas_image_height(image),jas_data);

		// transfer data

		k := ndata - image.height*image.width

		// #pragma omp parallel for private(ii,j)
		for ii := uzero; ii < image.height; ii++ {
			for j := uzero; j < image.width; j++ {
				// data[k++] = (((jas_data->rows_[ii][j])*bin_scale)+reference)*dec_scale;
				// data[k+j+ii*image.width] = (((jas_data->rows_[ii][j])*bin_scale)+reference)*dec_scale;
				data[k+j+ii*image.width] = 1.
			}
		}

		if bitmap_flag == 0 || bitmap_flag == 254 {
			k = ndata - image.height*image.width
			mask_pointer = sec[6][6:]
			mask = 0
			mask_pointer_index := 0
			for ii := uzero; ii < ndata; ii++ {
				if (ii & 7) == 0 {
					mask = mask_pointer[mask_pointer_index]
					mask_pointer_index++
				}
				// data[ii] = (mask & 128) ? data[k++] : UNDEFINED;
				if (mask & 128) != 0 {
					data[ii] = data[k]
					k++
				} else {
					data[ii] = UNDEFINED
				}
				mask <<= 1
			}
		} else if bitmap_flag != 255 {
			fatal_error_i("unknown bitmap: %d", bitmap_flag)
		}
		// jas_matrix_destroy(jas_data);
		// jas_stream_close(jpcstream);
		// jas_image_destroy(image);
		return nil
	} else if packing == 200 { // run length
		return fatal_error("unpk_run_length is not supported")
		// TODO: unpk_run_length
		// return unpk_run_length(sec, data, ndata)
	}
	return fatal_error_i("packing type %d not supported", packing)
}
