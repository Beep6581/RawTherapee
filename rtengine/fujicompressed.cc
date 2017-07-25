/* -*- C++ -*-
 * File: fujicompressed.cpp
 * Copyright (C) 2016 Alexey Danilchenko
 *
 * Adopted to LibRaw by Alex Tutubalin, lexa@lexa.ru
 * LibRaw Fujifilm/compressed decoder

LibRaw is free software; you can redistribute it and/or modify
it under the terms of the one of three licenses as you choose:

1. GNU LESSER GENERAL PUBLIC LICENSE version 2.1
   (See file LICENSE.LGPL provided in LibRaw distribution archive for details).

2. COMMON DEVELOPMENT AND DISTRIBUTION LICENSE (CDDL) Version 1.0
   (See file LICENSE.CDDL provided in LibRaw distribution archive for details).

 * Adopted to RawTherapee by Ingo Weyrich

 */

namespace {

int bitDiff (int value1, int value2)
{
    int decBits = 0;

    if ( value2 < value1 )
        while (decBits <= 12 && (value2 << ++decBits) < value1)
            ;

    return decBits;
}

}

void CLASS init_fuji_compr (struct fuji_compressed_params* info)
{
    int cur_val;
    char *qt;

    if ((fuji_block_width % 3 && fuji_raw_type == 16) || (fuji_block_width & 1 && fuji_raw_type == 0)) {
        derror();
    }

    info->q_table = (char *) malloc (32768);
    merror (info->q_table, "init_fuji_compr()");

    if (fuji_raw_type == 16) {
        info->line_width = (fuji_block_width * 2) / 3;
    } else {
        info->line_width = fuji_block_width >> 1;
    }

    info->q_point[0] = 0;
    info->q_point[1] = 0x12;
    info->q_point[2] = 0x43;
    info->q_point[3] = 0x114;
    info->q_point[4] = (1 << fuji_bits) - 1;
    info->min_value = 0x40;

    cur_val = -info->q_point[4];

    for (qt = info->q_table; cur_val <= info->q_point[4]; ++qt, ++cur_val) {
        if (cur_val <= -info->q_point[3]) {
            *qt = -4;
        } else if (cur_val <= -info->q_point[2]) {
            *qt = -3;
        } else if (cur_val <= -info->q_point[1]) {
            *qt = -2;
        } else if (cur_val < 0) {
            *qt = -1;
        } else if (cur_val == 0) {
            *qt = 0;
        } else if (cur_val < info->q_point[1]) {
            *qt = 1;
        } else if (cur_val < info->q_point[2]) {
            *qt = 2;
        } else if (cur_val < info->q_point[3]) {
            *qt = 3;
        } else {
            *qt = 4;
        }
    }

    // populting gradients
    if (info->q_point[4] == 0x3FFF) {
        info->total_values = 0x4000;
        info->raw_bits = 14;
        info->max_bits = 56;
        info->maxDiff = 256;
    } else if (info->q_point[4] == 0xFFF) {
        info->total_values = 4096;
        info->raw_bits = 12;
        info->max_bits = 48;
        info->maxDiff = 64;
    } else {
        derror();
    }
}

#define FUJI_BUF_SIZE 0x10000u

void CLASS fuji_fill_buffer (struct fuji_compressed_block *info)
{
    if (info->cur_pos >= info->cur_buf_size) {
        info->cur_pos = 0;
        info->cur_buf_offset += info->cur_buf_size;
#ifdef MYFILE_MMAP
        info->cur_buf_size = info->max_read_size;
        info->cur_buf = fdata(info->cur_buf_offset, info->input);
#else
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            fseek (info->input, info->cur_buf_offset, SEEK_SET);
            info->cur_buf_size = fread (info->cur_buf, 1, std::min (info->max_read_size, FUJI_BUF_SIZE), info->input);
        }
#endif
        if (info->cur_buf_size < 1) { // nothing read
            if (info->fillbytes > 0) {
                int ls = std::max (1, std::min (info->fillbytes, (int)FUJI_BUF_SIZE));
                memset (info->cur_buf, 0, ls);
                info->fillbytes -= ls;
            }
        }

        info->max_read_size -= info->cur_buf_size;
    }
}

void CLASS init_fuji_block (struct fuji_compressed_block* info, const struct fuji_compressed_params *params, INT64 raw_offset, unsigned dsize)
{
    info->linealloc = (ushort*)calloc (sizeof (ushort), _ltotal * (params->line_width + 2));
    merror (info->linealloc, "init_fuji_block()");

    info->input = ifp;
    INT64 fsize = info->input->size;
    info->max_read_size = std::min (unsigned (fsize - raw_offset), dsize + 16); // Data size may be incorrect?
    info->fillbytes = 1;

    info->linebuf[_R0] = info->linealloc;

    for (int i = _R1; i <= _B4; i++) {
        info->linebuf[i] = info->linebuf[i - 1] + params->line_width + 2;
    }

    // init buffer
#ifndef MYFILE_MMAP
    info->cur_buf = (uchar*)malloc (FUJI_BUF_SIZE);
    merror (info->cur_buf, "init_fuji_block()");
#endif
    info->cur_bit = 0;
    info->cur_pos = 0;
    info->cur_buf_offset = raw_offset;

    for (int j = 0; j < 3; j++)
        for (int i = 0; i < 41; i++) {
            info->grad_even[j][i].value1 = params->maxDiff;
            info->grad_even[j][i].value2 = 1;
            info->grad_odd[j][i].value1 = params->maxDiff;
            info->grad_odd[j][i].value2 = 1;
        }

    info->cur_buf_size = 0;
    fuji_fill_buffer (info);
}

void CLASS copy_line_to_xtrans (struct fuji_compressed_block* info, int cur_line, int cur_block, int cur_block_width)
{
    ushort *lineBufB[3];
    ushort *lineBufG[6];
    ushort *lineBufR[3];
    ushort* line_buf;
    int index;

    int offset = fuji_block_width * cur_block + 6 * raw_width * cur_line;
    ushort* raw_block_data = raw_image + offset;
    int row_count = 0;

    for (int i = 0; i < 3; i++) {
        lineBufR[i] = info->linebuf[_R2 + i] + 1;
        lineBufB[i] = info->linebuf[_B2 + i] + 1;
    }

    for (int i = 0; i < 6; i++) {
        lineBufG[i] = info->linebuf[_G2 + i] + 1;
    }

    while (row_count < 6) {
        unsigned pixel_count = 0;

        while (static_cast<int>(pixel_count) < cur_block_width) {
            switch (xtrans_abs[row_count][ (pixel_count % 6)]) {
                case 1:     // green
                    line_buf = lineBufG[row_count];
                    break;

                case 0:     // red
                    line_buf = lineBufR[row_count >> 1];
                    break;

                case 2:     // blue
                default:
                    line_buf = lineBufB[row_count >> 1];
                    break;
            }

            index = (((pixel_count * 2 / 3) & 0x7FFFFFFE) | ((pixel_count % 3) & 1)) + ((pixel_count % 3) >> 1);
            raw_block_data[pixel_count] = line_buf[index];

            ++pixel_count;
        }

        ++row_count;
        raw_block_data += raw_width;
    }
}

void CLASS copy_line_to_bayer (struct fuji_compressed_block *info, int cur_line, int cur_block, int cur_block_width)
{
    ushort *lineBufB[3];
    ushort *lineBufG[6];
    ushort *lineBufR[3];
    ushort *line_buf;

    int fuji_bayer[2][2];

    for (int r = 0; r < 2; r++)
        for (int c = 0; c < 2; c++) {
            fuji_bayer[r][c] = FC (r, c);    // We'll downgrade G2 to G below
        }

    int offset = fuji_block_width * cur_block + 6 * raw_width * cur_line;
    ushort *raw_block_data = raw_image + offset;
    int row_count = 0;

    for (int i = 0; i < 3; i++) {
        lineBufR[i] = info->linebuf[_R2 + i] + 1;
        lineBufB[i] = info->linebuf[_B2 + i] + 1;
    }

    for (int i = 0; i < 6; i++) {
        lineBufG[i] = info->linebuf[_G2 + i] + 1;
    }

    while (row_count < 6) {
        unsigned pixel_count = 0;

        while (static_cast<int>(pixel_count) < cur_block_width) {
            switch (fuji_bayer[row_count & 1][pixel_count & 1]) {
                case 0: // red
                    line_buf = lineBufR[row_count >> 1];
                    break;

                case 1:  // green
                case 3:  // second green
                default: // to make static analyzer happy
                    line_buf = lineBufG[row_count];
                    break;

                case 2: // blue
                    line_buf = lineBufB[row_count >> 1];
                    break;
            }

            raw_block_data[pixel_count] = line_buf[pixel_count >> 1];
            ++pixel_count;
        }

        ++row_count;
        raw_block_data += raw_width;
    }
}


#define fuji_quant_gradient(i,v1,v2) (9*i->q_table[i->q_point[4]+(v1)] + i->q_table[i->q_point[4]+(v2)])

inline void CLASS fuji_zerobits (struct fuji_compressed_block* info, int *count)
{
    uchar zero = 0;
    *count = 0;

    while (zero == 0) {
        zero = (info->cur_buf[info->cur_pos] >> (7 - info->cur_bit)) & 1;
        info->cur_bit++;
        info->cur_bit &= 7;

        if (!info->cur_bit) {
            ++info->cur_pos;
#ifndef MYFILE_MMAP
            fuji_fill_buffer (info);
#endif
        }

        if (zero) {
            break;
        }

        ++*count;
    }
}

inline void CLASS fuji_read_code (struct fuji_compressed_block* info, int *data, int bits_to_read)
{
    *data = 0;

    if (!bits_to_read) {
        return;
    }

    uchar bits_left = bits_to_read;
    uchar bits_left_in_byte = 8 - (info->cur_bit & 7);

    if (bits_to_read >= bits_left_in_byte) {
        do {
            *data <<= bits_left_in_byte;
            bits_left -= bits_left_in_byte;
            *data |= info->cur_buf[info->cur_pos] & ((1 << bits_left_in_byte) - 1);
            ++info->cur_pos;
#ifndef MYFILE_MMAP
            fuji_fill_buffer (info);
#endif
            bits_left_in_byte = 8;
        } while (bits_left >= 8);
    }

    if (!bits_left) {
        info->cur_bit = (8 - (bits_left_in_byte & 7)) & 7;
        return;
    }

    *data <<= bits_left;
    bits_left_in_byte -= bits_left;
    *data |= ((1 << bits_left) - 1) & ((unsigned)info->cur_buf[info->cur_pos] >> bits_left_in_byte);
    info->cur_bit = (8 - (bits_left_in_byte & 7)) & 7;
}

int CLASS fuji_decode_sample_even (struct fuji_compressed_block* info, const struct fuji_compressed_params * params, ushort* line_buf, int pos, struct int_pair* grads)
{
    int interp_val = 0;
    int errcnt = 0;

    int sample = 0, code = 0;
    ushort* line_buf_cur = line_buf + pos;
    int Rb = line_buf_cur[-2 - params->line_width];
    int Rc = line_buf_cur[-3 - params->line_width];
    int Rd = line_buf_cur[-1 - params->line_width];
    int Rf = line_buf_cur[-4 - 2 * params->line_width];

    int grad, gradient, diffRcRb, diffRfRb, diffRdRb;

    grad = fuji_quant_gradient (params, Rb - Rf, Rc - Rb);
    gradient = std::abs (grad);
    diffRcRb = std::abs (Rc - Rb);
    diffRfRb = std::abs (Rf - Rb);
    diffRdRb = std::abs (Rd - Rb);

    if ( diffRcRb > diffRfRb && diffRcRb > diffRdRb ) {
        interp_val = Rf + Rd + 2 * Rb;
    } else if ( diffRdRb > diffRcRb && diffRdRb > diffRfRb ) {
        interp_val = Rf + Rc + 2 * Rb;
    } else {
        interp_val = Rd + Rc + 2 * Rb;
    }


    fuji_zerobits (info, &sample);

    if (sample < params->max_bits - params->raw_bits - 1) {
        int decBits = bitDiff (grads[gradient].value1, grads[gradient].value2);
        fuji_read_code (info, &code, decBits);
        code += sample << decBits;
    } else {
        fuji_read_code (info, &code, params->raw_bits);
        code++;
    }

    if (code < 0 || code >= params->total_values) {
        errcnt++;
    }

    if (code & 1) {
        code = -1 - code / 2;
    } else {
        code /= 2;
    }

    grads[gradient].value1 += std::abs (code);

    if (grads[gradient].value2 == params->min_value ) {
        grads[gradient].value1 >>= 1;
        grads[gradient].value2 >>= 1;
    }

    grads[gradient].value2++;

    if (grad < 0) {
        interp_val = (interp_val >> 2) - code;
    } else {
        interp_val = (interp_val >> 2) + code;
    }

    if ( interp_val < 0 ) {
        interp_val += params->total_values;
    } else if (interp_val > params->q_point[4]) {
        interp_val -= params->total_values;
    }

    if ( interp_val >= 0 ) {
        line_buf_cur[0] = std::min (interp_val, params->q_point[4]);
    } else {
        line_buf_cur[0] = 0;
    }

    return errcnt;
}

int CLASS fuji_decode_sample_odd (struct fuji_compressed_block* info, const struct fuji_compressed_params * params, ushort* line_buf, int pos, struct int_pair* grads)
{
    int interp_val = 0;
    int errcnt = 0;

    int sample = 0, code = 0;
    ushort* line_buf_cur = line_buf + pos;
    int Ra = line_buf_cur[-1];
    int Rb = line_buf_cur[-2 - params->line_width];
    int Rc = line_buf_cur[-3 - params->line_width];
    int Rd = line_buf_cur[-1 - params->line_width];
    int Rg = line_buf_cur[1];

    int grad, gradient;

    grad = fuji_quant_gradient (params, Rb - Rc, Rc - Ra);
    gradient = std::abs (grad);

    if ((Rb > Rc && Rb > Rd) || (Rb < Rc && Rb < Rd)) {
        interp_val = (Rg + Ra + 2 * Rb) >> 2;
    } else {
        interp_val = (Ra + Rg) >> 1;
    }

    fuji_zerobits (info, &sample);

    if (sample < params->max_bits - params->raw_bits - 1) {
        int decBits = bitDiff (grads[gradient].value1, grads[gradient].value2);
        fuji_read_code (info, &code, decBits);
        code += sample << decBits;
    } else {
        fuji_read_code (info, &code, params->raw_bits);
        code++;
    }

    if (code < 0 || code >= params->total_values) {
        errcnt++;
    }

    if (code & 1) {
        code = -1 - code / 2;
    } else {
        code /= 2;
    }

    grads[gradient].value1 += std::abs (code);

    if (grads[gradient].value2 == params->min_value) {
        grads[gradient].value1 >>= 1;
        grads[gradient].value2 >>= 1;
    }

    grads[gradient].value2++;

    if (grad < 0) {
        interp_val -= code;
    } else {
        interp_val += code;
    }

    if ( interp_val < 0 ) {
        interp_val += params->total_values;
    } else if (interp_val > params->q_point[4]) {
        interp_val -= params->total_values;
    }

    if ( interp_val >= 0 ) {
        line_buf_cur[0] = std::min (interp_val, params->q_point[4]);
    } else {
        line_buf_cur[0] = 0;
    }

    return errcnt;
}

void CLASS fuji_decode_interpolation_even (int line_width, ushort* line_buf, int pos)
{
    ushort* line_buf_cur = line_buf + pos;
    int Rb = line_buf_cur[-2 - line_width];
    int Rc = line_buf_cur[-3 - line_width];
    int Rd = line_buf_cur[-1 - line_width];
    int Rf = line_buf_cur[-4 - 2 * line_width];
    int diffRcRb = std::abs (Rc - Rb);
    int diffRfRb = std::abs (Rf - Rb);
    int diffRdRb = std::abs (Rd - Rb);

    if ( diffRcRb > diffRfRb && diffRcRb > diffRdRb ) {
        *line_buf_cur = (Rf + Rd + 2 * Rb) >> 2;
    } else if ( diffRdRb > diffRcRb && diffRdRb > diffRfRb ) {
        *line_buf_cur = (Rf + Rc + 2 * Rb) >> 2;
    } else {
        *line_buf_cur = (Rd + Rc + 2 * Rb) >> 2;
    }
}

void CLASS fuji_extend_generic (ushort *linebuf[_ltotal], int line_width, int start, int end)
{
    for (int i = start; i <= end; i++) {
        linebuf[i][0]             = linebuf[i - 1][1];
        linebuf[i][line_width + 1] = linebuf[i - 1][line_width];
    }
}

void CLASS fuji_extend_red (ushort *linebuf[_ltotal], int line_width)
{
    fuji_extend_generic (linebuf, line_width, _R2, _R4);
}

void CLASS fuji_extend_green (ushort *linebuf[_ltotal], int line_width)
{
    fuji_extend_generic (linebuf, line_width, _G2, _G7);
}

void CLASS fuji_extend_blue (ushort *linebuf[_ltotal], int line_width)
{
    fuji_extend_generic (linebuf, line_width, _B2, _B4);
}

void CLASS xtrans_decode_block (struct fuji_compressed_block* info, const struct fuji_compressed_params *params, int cur_line)
{
    int r_even_pos = 0, r_odd_pos = 1;
    int g_even_pos = 0, g_odd_pos = 1;
    int b_even_pos = 0, b_odd_pos = 1;

    int errcnt = 0;

    const int line_width = params->line_width;

    while (g_even_pos < line_width || g_odd_pos < line_width) {
        if (g_even_pos < line_width) {
            fuji_decode_interpolation_even (line_width, info->linebuf[_R2] + 1, r_even_pos);
            r_even_pos += 2;
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_G2] + 1, g_even_pos, info->grad_even[0]);
            g_even_pos += 2;
        }

        if (g_even_pos > 8) {
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_R2] + 1, r_odd_pos, info->grad_odd[0]);
            r_odd_pos += 2;
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_G2] + 1, g_odd_pos, info->grad_odd[0]);
            g_odd_pos += 2;
        }
    }

    fuji_extend_red (info->linebuf, line_width);
    fuji_extend_green (info->linebuf, line_width);

    g_even_pos = 0, g_odd_pos = 1;

    while (g_even_pos < line_width || g_odd_pos < line_width) {
        if (g_even_pos < line_width) {
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_G3] + 1, g_even_pos, info->grad_even[1]);
            g_even_pos += 2;
            fuji_decode_interpolation_even (line_width, info->linebuf[_B2] + 1, b_even_pos);
            b_even_pos += 2;
        }

        if (g_even_pos > 8) {
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_G3] + 1, g_odd_pos, info->grad_odd[1]);
            g_odd_pos += 2;
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_B2] + 1, b_odd_pos, info->grad_odd[1]);
            b_odd_pos += 2;
        }
    }

    fuji_extend_green (info->linebuf, line_width);
    fuji_extend_blue (info->linebuf, line_width);

    r_even_pos = 0, r_odd_pos = 1;
    g_even_pos = 0, g_odd_pos = 1;

    while (g_even_pos < line_width || g_odd_pos < line_width) {
        if (g_even_pos < line_width) {
            if (r_even_pos & 3) {
                errcnt += fuji_decode_sample_even (info, params, info->linebuf[_R3] + 1, r_even_pos, info->grad_even[2]);
            } else {
                fuji_decode_interpolation_even (line_width, info->linebuf[_R3] + 1, r_even_pos);
            }

            r_even_pos += 2;
            fuji_decode_interpolation_even (line_width, info->linebuf[_G4] + 1, g_even_pos);
            g_even_pos += 2;
        }

        if (g_even_pos > 8) {
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_R3] + 1, r_odd_pos, info->grad_odd[2]);
            r_odd_pos += 2;
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_G4] + 1, g_odd_pos, info->grad_odd[2]);
            g_odd_pos += 2;
        }
    }

    fuji_extend_red (info->linebuf, line_width);
    fuji_extend_green (info->linebuf, line_width);

    g_even_pos = 0, g_odd_pos = 1;
    b_even_pos = 0, b_odd_pos = 1;

    while (g_even_pos < line_width || g_odd_pos < line_width) {
        if (g_even_pos < line_width) {
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_G5] + 1, g_even_pos, info->grad_even[0]);
            g_even_pos += 2;

            if ((b_even_pos & 3) == 2) {
                fuji_decode_interpolation_even (line_width, info->linebuf[_B3] + 1, b_even_pos);
            } else {
                errcnt += fuji_decode_sample_even (info, params, info->linebuf[_B3] + 1, b_even_pos, info->grad_even[0]);
            }

            b_even_pos += 2;
        }

        if (g_even_pos > 8) {
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_G5] + 1, g_odd_pos, info->grad_odd[0]);
            g_odd_pos += 2;
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_B3] + 1, b_odd_pos, info->grad_odd[0]);
            b_odd_pos += 2;
        }
    }

    fuji_extend_green (info->linebuf, line_width);
    fuji_extend_blue (info->linebuf, line_width);

    r_even_pos = 0, r_odd_pos = 1;
    g_even_pos = 0, g_odd_pos = 1;

    while (g_even_pos < line_width || g_odd_pos < line_width) {
        if (g_even_pos < line_width) {
            if ((r_even_pos & 3) == 2) {
                fuji_decode_interpolation_even (line_width, info->linebuf[_R4] + 1, r_even_pos);
            } else {
                errcnt += fuji_decode_sample_even (info, params, info->linebuf[_R4] + 1, r_even_pos, info->grad_even[1]);
            }

            r_even_pos += 2;
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_G6] + 1, g_even_pos, info->grad_even[1]);
            g_even_pos += 2;
        }

        if (g_even_pos > 8) {
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_R4] + 1, r_odd_pos, info->grad_odd[1]);
            r_odd_pos += 2;
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_G6] + 1, g_odd_pos, info->grad_odd[1]);
            g_odd_pos += 2;
        }
    }

    fuji_extend_red (info->linebuf, line_width);
    fuji_extend_green (info->linebuf, line_width);

    g_even_pos = 0, g_odd_pos = 1;
    b_even_pos = 0, b_odd_pos = 1;

    while (g_even_pos < line_width || g_odd_pos < line_width) {
        if (g_even_pos < line_width) {
            fuji_decode_interpolation_even (line_width, info->linebuf[_G7] + 1, g_even_pos);
            g_even_pos += 2;

            if (b_even_pos & 3) {
                errcnt += fuji_decode_sample_even (info, params, info->linebuf[_B4] + 1, b_even_pos, info->grad_even[2]);
            } else {
                fuji_decode_interpolation_even (line_width, info->linebuf[_B4] + 1, b_even_pos);
            }

            b_even_pos += 2;
        }

        if (g_even_pos > 8) {
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_G7] + 1, g_odd_pos, info->grad_odd[2]);
            g_odd_pos += 2;
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_B4] + 1, b_odd_pos, info->grad_odd[2]);
            b_odd_pos += 2;
        }
    }

    fuji_extend_green (info->linebuf, line_width);
    fuji_extend_blue (info->linebuf, line_width);

    if (errcnt) {
        derror();
    }
}

void CLASS fuji_bayer_decode_block (struct fuji_compressed_block *info, const struct fuji_compressed_params *params,
                                    int cur_line)
{
    int r_even_pos = 0, r_odd_pos = 1;
    int g_even_pos = 0, g_odd_pos = 1;
    int b_even_pos = 0, b_odd_pos = 1;

    int errcnt = 0;

    const int line_width = params->line_width;

    while (g_even_pos < line_width || g_odd_pos < line_width) {
        if (g_even_pos < line_width) {
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_R2] + 1, r_even_pos, info->grad_even[0]);
            r_even_pos += 2;
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_G2] + 1, g_even_pos, info->grad_even[0]);
            g_even_pos += 2;
        }

        if (g_even_pos > 8) {
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_R2] + 1, r_odd_pos, info->grad_odd[0]);
            r_odd_pos += 2;
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_G2] + 1, g_odd_pos, info->grad_odd[0]);
            g_odd_pos += 2;
        }
    }

    fuji_extend_red (info->linebuf, line_width);
    fuji_extend_green (info->linebuf, line_width);

    g_even_pos = 0, g_odd_pos = 1;

    while (g_even_pos < line_width || g_odd_pos < line_width) {
        if (g_even_pos < line_width) {
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_G3] + 1, g_even_pos, info->grad_even[1]);
            g_even_pos += 2;
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_B2] + 1, b_even_pos, info->grad_even[1]);
            b_even_pos += 2;
        }

        if (g_even_pos > 8) {
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_G3] + 1, g_odd_pos, info->grad_odd[1]);
            g_odd_pos += 2;
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_B2] + 1, b_odd_pos, info->grad_odd[1]);
            b_odd_pos += 2;
        }
    }

    fuji_extend_green (info->linebuf, line_width);
    fuji_extend_blue (info->linebuf, line_width);

    r_even_pos = 0, r_odd_pos = 1;
    g_even_pos = 0, g_odd_pos = 1;

    while (g_even_pos < line_width || g_odd_pos < line_width) {
        if (g_even_pos < line_width) {
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_R3] + 1, r_even_pos, info->grad_even[2]);
            r_even_pos += 2;
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_G4] + 1, g_even_pos, info->grad_even[2]);
            g_even_pos += 2;
        }

        if (g_even_pos > 8) {
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_R3] + 1, r_odd_pos, info->grad_odd[2]);
            r_odd_pos += 2;
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_G4] + 1, g_odd_pos, info->grad_odd[2]);
            g_odd_pos += 2;
        }
    }

    fuji_extend_red (info->linebuf, line_width);
    fuji_extend_green (info->linebuf, line_width);

    g_even_pos = 0, g_odd_pos = 1;
    b_even_pos = 0, b_odd_pos = 1;

    while (g_even_pos < line_width || g_odd_pos < line_width) {
        if (g_even_pos < line_width) {
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_G5] + 1, g_even_pos, info->grad_even[0]);
            g_even_pos += 2;
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_B3] + 1, b_even_pos, info->grad_even[0]);
            b_even_pos += 2;
        }

        if (g_even_pos > 8) {
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_G5] + 1, g_odd_pos, info->grad_odd[0]);
            g_odd_pos += 2;
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_B3] + 1, b_odd_pos, info->grad_odd[0]);
            b_odd_pos += 2;
        }
    }

    fuji_extend_green (info->linebuf, line_width);
    fuji_extend_blue (info->linebuf, line_width);

    r_even_pos = 0, r_odd_pos = 1;
    g_even_pos = 0, g_odd_pos = 1;

    while (g_even_pos < line_width || g_odd_pos < line_width) {
        if (g_even_pos < line_width) {
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_R4] + 1, r_even_pos, info->grad_even[1]);
            r_even_pos += 2;
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_G6] + 1, g_even_pos, info->grad_even[1]);
            g_even_pos += 2;
        }

        if (g_even_pos > 8) {
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_R4] + 1, r_odd_pos, info->grad_odd[1]);
            r_odd_pos += 2;
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_G6] + 1, g_odd_pos, info->grad_odd[1]);
            g_odd_pos += 2;
        }
    }

    fuji_extend_red (info->linebuf, line_width);
    fuji_extend_green (info->linebuf, line_width);

    g_even_pos = 0, g_odd_pos = 1;
    b_even_pos = 0, b_odd_pos = 1;

    while (g_even_pos < line_width || g_odd_pos < line_width) {
        if (g_even_pos < line_width) {
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_G7] + 1, g_even_pos, info->grad_even[2]);
            g_even_pos += 2;
            errcnt += fuji_decode_sample_even (info, params, info->linebuf[_B4] + 1, b_even_pos, info->grad_even[2]);
            b_even_pos += 2;
        }

        if (g_even_pos > 8) {
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_G7] + 1, g_odd_pos, info->grad_odd[2]);
            g_odd_pos += 2;
            errcnt += fuji_decode_sample_odd (info, params, info->linebuf[_B4] + 1, b_odd_pos, info->grad_odd[2]);
            b_odd_pos += 2;
        }
    }

    fuji_extend_green (info->linebuf, line_width);
    fuji_extend_blue (info->linebuf, line_width);

    if (errcnt) {
        derror();
    }
}

void CLASS fuji_decode_strip (const struct fuji_compressed_params* info_common, int cur_block, INT64 raw_offset, unsigned dsize)
{
    int cur_block_width, cur_line;
    unsigned line_size;
    struct fuji_compressed_block info;

    init_fuji_block (&info, info_common, raw_offset, dsize);
    line_size = sizeof (ushort) * (info_common->line_width + 2);

    cur_block_width = fuji_block_width;

    if (cur_block + 1 == fuji_total_blocks) {
        cur_block_width = raw_width - cur_block * fuji_block_width;
    }

    struct i_pair {
        int a, b;
    };

    const i_pair mtable[6] = { {_R0, _R3}, {_R1, _R4}, {_G0, _G6}, {_G1, _G7}, {_B0, _B3}, {_B1, _B4}},
    ztable[3] = {{_R2, 3}, {_G2, 6}, {_B2, 3}};

    for  (cur_line = 0; cur_line < fuji_total_lines; cur_line++) {
        if (fuji_raw_type == 16) {
            xtrans_decode_block (&info, info_common, cur_line);
        } else {
            fuji_bayer_decode_block (&info, info_common, cur_line);
        }

        // copy data from line buffers and advance
        for (int i = 0; i < 6; i++) {
            memcpy (info.linebuf[mtable[i].a], info.linebuf[mtable[i].b], line_size);
        }

        if (fuji_raw_type == 16) {
            copy_line_to_xtrans (&info, cur_line, cur_block, cur_block_width);
        } else {
            copy_line_to_bayer (&info, cur_line, cur_block, cur_block_width);
        }

        for (int i = 0; i < 3; i++) {
            memset (info.linebuf[ztable[i].a], 0, ztable[i].b * line_size);
            info.linebuf[ztable[i].a][0]                    = info.linebuf[ztable[i].a - 1][1];
            info.linebuf[ztable[i].a][info_common->line_width + 1] = info.linebuf[ztable[i].a - 1][info_common->line_width];
        }
    }

    // release data
    free (info.linealloc);
#ifndef MYFILE_MMAP
    free (info.cur_buf);
#endif
}

static unsigned sgetn (int n, uchar *s)
{
    unsigned result = 0;

    while (n-- > 0) {
        result = (result << 8) | (*s++);
    }

    return result;
}

void CLASS fuji_compressed_load_raw()
{
    struct fuji_compressed_params common_info;
    int cur_block;
    unsigned *block_sizes;
    INT64 raw_offset, *raw_block_offsets;
    //struct fuji_compressed_block info;

    init_fuji_compr (&common_info);

    // read block sizes
    block_sizes = (unsigned*) malloc (sizeof (unsigned) * fuji_total_blocks);
    merror (block_sizes, "fuji_compressed_load_raw()");
    raw_block_offsets = (INT64*) malloc (sizeof (INT64) * fuji_total_blocks);
    merror (raw_block_offsets, "fuji_compressed_load_raw()");

    raw_offset = sizeof (unsigned) * fuji_total_blocks;

    if (raw_offset & 0xC) {
        raw_offset +=  0x10 - (raw_offset & 0xC);
    }

    raw_offset += data_offset;

    fseek (ifp, data_offset, SEEK_SET);
    fread (block_sizes, 1, sizeof (unsigned)*fuji_total_blocks, ifp);

    raw_block_offsets[0] = raw_offset;

    // calculating raw block offsets
    for (cur_block = 0; cur_block < fuji_total_blocks; cur_block++) {
        unsigned bsize = sgetn (4, (uchar *) (block_sizes + cur_block));
        block_sizes[cur_block] = bsize;
    }

    for (cur_block = 1; cur_block < fuji_total_blocks; cur_block++) {
        raw_block_offsets[cur_block] = raw_block_offsets[cur_block - 1] + block_sizes[cur_block - 1] ;
    }

    fuji_decode_loop (&common_info, fuji_total_blocks, raw_block_offsets, block_sizes);

    free (block_sizes);
    free (raw_block_offsets);
    free (common_info.q_table);
}

void CLASS fuji_decode_loop (const struct fuji_compressed_params* common_info, int count, INT64* raw_block_offsets, unsigned *block_sizes)
{

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1) // dynamic scheduling is faster if count > number of cores (e.g. count for GFX 50S is 12)
#endif

    for (int cur_block = 0; cur_block < count ; cur_block++) {
        fuji_decode_strip (common_info, cur_block, raw_block_offsets[cur_block], block_sizes[cur_block]);
    }
}


void CLASS parse_fuji_compressed_header()
{

    uchar header[16];

    ushort signature;
    uchar  version;
    uchar  h_raw_type;
    uchar  h_raw_bits;
    ushort h_raw_height;
    ushort h_raw_rounded_width;
    ushort h_raw_width;
    ushort h_block_size;
    uchar  h_blocks_in_row;
    ushort h_total_lines;

    fseek (ifp, data_offset, SEEK_SET);
    fread (header, 1, sizeof (header), ifp);

    signature = sgetn (2, header);
    version = header[2];
    h_raw_type = header[3];
    h_raw_bits = header[4];
    h_raw_height = sgetn (2, header + 5);
    h_raw_rounded_width = sgetn (2, header + 7);
    h_raw_width = sgetn (2, header + 9);
    h_block_size = sgetn (2, header + 11);
    h_blocks_in_row = header[13];
    h_total_lines = sgetn (2, header + 14);


    // general validation
    if (signature != 0x4953
            || version != 1
            || h_raw_height > 0x3000
            || h_raw_height < 6
            || h_raw_height % 6
            || h_raw_width > 0x3000
            || h_raw_width < 0x300
            || h_raw_width % 24
            || h_raw_rounded_width > 0x3000
            || h_block_size != 0x300
            || h_raw_rounded_width < h_block_size
            || h_raw_rounded_width % h_block_size
            || h_raw_rounded_width - h_raw_width >= h_block_size
            || h_blocks_in_row > 0x10
            || h_blocks_in_row == 0
            || h_blocks_in_row != h_raw_rounded_width / h_block_size
            || h_total_lines > 0x800
            || h_total_lines == 0
            || h_total_lines != h_raw_height / 6
            || (h_raw_bits != 12 && h_raw_bits != 14)
            || (h_raw_type != 16 && h_raw_type != 0)) {
        xtransCompressed = false;
        return;
    }

    // modify data
    fuji_total_lines  = h_total_lines;
    fuji_total_blocks = h_blocks_in_row;
    fuji_block_width  = h_block_size;
    fuji_bits         = h_raw_bits;
    fuji_raw_type     = h_raw_type;
    raw_width         = h_raw_width;
    raw_height        = h_raw_height;
    data_offset += 16;
    load_raw = &CLASS fuji_compressed_load_raw;
}
