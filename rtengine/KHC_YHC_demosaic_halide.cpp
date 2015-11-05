//
// Created by Lucas Cooper on 3/11/2015.
//

#include "KHC_YHC_demosaic_halide.h"

using namespace Halide;

Func make_demosaic_func(ImageParam param, Type out_type) {
    Var x, y, c;
    x = Var("x");
    y = Var("y");
    c = Var("c");
    Expr alpha = 1.0f;

    Func input = Func("input");
    input(x, y) = cast(Float(64), param(
            clamp(x, 0, param.width()-1),
            clamp(y, 0, param.height()-1)
    ));

    //#### Horizontal ####
    Func cd_horizontal = Func("cd_horizontal");
    cd_horizontal(x,y) = input(x,y) - (input(x-1,y) + input(x+1,y))/2;

    Func cd_gradient_horizontal_raw = Func("cd_gradient_horizontal_raw");
    cd_gradient_horizontal_raw(x,y) = (
                                              absd(cd_horizontal(x,y), cd_horizontal(x+1,y)) +
                                              absd(cd_horizontal(x,y), cd_horizontal(x-1,y)))/2.0f;

    Func cd_gradient_horizontal = Func("cd_gradient_horizontal");
    cd_gradient_horizontal(x,y) = (cd_gradient_horizontal_raw(x-1,y) + cd_gradient_horizontal_raw(x,y) + cd_gradient_horizontal_raw(x+1,y))/3.0f;

    Func ci_gradient_plane_horizontal = Func("ci_gradient_plane_horizontal");
    ci_gradient_plane_horizontal(x,y) = absd(input(x, y), input(x+2, y));

    Func cd_gradient_plane_horizontal = Func("cd_gradient_plane_horizontal");
    cd_gradient_plane_horizontal(x,y) = cd_gradient_horizontal(x+1, y);

    Func ig_horizontal = Func("ig_horizontal");
    ig_horizontal(x,y) = ci_gradient_plane_horizontal(x,y) + alpha * (
            2 * cd_gradient_plane_horizontal(x,y) +
            cd_gradient_plane_horizontal(x,y-1) +
            cd_gradient_plane_horizontal(x,y+1));

    //#### Vertical ####

    Func cd_vertical = Func("cd_vertical");
    cd_vertical(x,y) = input(x,y) - (input(x,y-1) + input(x,y+1))/2;

    Func cd_gradient_vertical_raw = Func("cd_gradient_vertical_raw");
    cd_gradient_vertical_raw(x,y) = (
                                            absd(cd_vertical(x,y), cd_vertical(x,y+1)) +
                                            absd(cd_vertical(x,y), cd_vertical(x,y-1)))/2.0f;

    Func cd_gradient_vertical = Func("cd_gradient_vertical");
    cd_gradient_vertical(x,y) = (cd_gradient_vertical_raw(x,y-1) + cd_gradient_vertical_raw(x,y) + cd_gradient_vertical_raw(x,y+1))/3.0f;

    Func ci_gradient_plane_vertical = Func("ci_gradient_plane_vertical");
    ci_gradient_plane_vertical(x,y) = absd(input(x, y), input(x, y+2));

    Func cd_gradient_plane_vertical = Func("cd_gradient_plane_vertical");
    cd_gradient_plane_vertical(x,y) = cd_gradient_vertical(x, y+1);


    Func ig_vertical = Func("ig_vertical");
    ig_vertical(x,y) = ci_gradient_plane_vertical(x,y) + alpha * (
            2 * cd_gradient_plane_vertical(x,y) +
            cd_gradient_plane_vertical(x-1,y) +
            cd_gradient_plane_vertical(x+1,y));

    //#### Interpolation ####
    Func green_horizontal = Func("green_horizontal");
    green_horizontal(x,y) = ((input(x-1,y) + input(x+1,y))/2.0f) + (2.0f * input(x,y) - input(x-2,y) - input(x+2, y))/4.0f;

    Func green_vertical = Func("green_vertical");
    green_vertical(x,y) = ((input(x,y-1) + input(x,y+1))/2.0f) + (2.0f * input(x,y) - input(x,y-2) - input(x, y+2))/4.0f;

    Func green_neither = Func("green_neither");
    green_neither(x,y) = (green_horizontal(x,y) + green_vertical(x,y))/2.0f;


    Func delta_h = Func("delta_h");
    delta_h(x,y) = ig_horizontal(x,y) + ig_horizontal(x-2,y);

    Func delta_v = Func("delta_v");
    delta_v(x,y) = ig_vertical(x,y) + ig_vertical(x,y-2);

    Func max_factor = Func("max_factor");
    max_factor(x,y) = max(delta_h(x,y)/delta_v(x,y), delta_v(x,y)/delta_h(x,y));

    float threshold = 1.7f;
    Expr L = 3.0f;


    Func go_to_pass1 = Func("go_to_pass1");
    go_to_pass1(x,y) = select(
            max_factor(x,y) == 1 or max_factor(x,y) > threshold, 1,
            0
    );


    Func pass1 = Func("pass1");
    pass1(x,y) = select(
            delta_h(x,y) == delta_v(x,y), green_neither(x,y),
            delta_h(x,y) > delta_v(x,y), green_vertical(x,y),
            green_horizontal(x,y));


    Func pmn_h = Func("pmn_h");
    pmn_h(x,y) = select(
            go_to_pass1(x,y) == 1, pass1(x,y) - input(x,y),
            green_horizontal(x,y) - input(x,y)
    );


    Func pmn_v = Func("pmn_v");
    pmn_v(x,y) = select(
            go_to_pass1(x,y) == 1, pass1(x,y) - input(x,y),
            green_vertical(x,y) - input(x,y)
    );


    Func pmn_d = Func("pmn_d");
    pmn_d(x,y) = select(
            go_to_pass1(x,y) == 1, pass1(x,y) - input(x,y),
            green_neither(x,y) - input(x,y)
    );

    RDom tDom = RDom(-L,L+1);

    Func phi_h = Func("phi_h");
    phi_h(x,y) = sum(absd(pmn_h(x,y), pmn_h(x-2 * tDom.x,y)));

    Func phi_v = Func("phi_v");
    phi_v(x,y) = sum(absd(pmn_v(x,y), pmn_v(x,y-2 * tDom.x)));

    Func phi_d = Func("phi_d");
    phi_d(x,y) = (
                         sum(absd(pmn_d(x,y), pmn_d(x-2 * tDom.x, y))) +
                         sum(absd(pmn_d(x,y), pmn_d(x, y-2 * tDom.x))))/2.0f;

    Func smallest_phi = Func("smallest_phi");
    smallest_phi(x,y) = Halide::min(phi_h(x,y), Halide::min(phi_v(x,y), phi_d(x,y)));

    Func pass2 = Func("pass2");
    pass2(x,y) = select(
            go_to_pass1(x,y) == 1, pass1(x,y),
            smallest_phi(x,y) == phi_h(x,y), green_horizontal(x,y),
            smallest_phi(x,y) == phi_v(x,y), green_vertical(x,y),
            green_neither(x,y)
    );


    Func green_unenhanced = Func("green_unenhanced");
    green_unenhanced(x,y) = select(
            x%2+y%2 == 1, input(x,y),
            pass2(x,y));

    Func w_e = Func("w_e");
    w_e(x,y) = 1/ig_horizontal(x,y);
    Func w_w = Func("w_w");
    w_w(x,y) = 1/ig_horizontal(x-2,y);
    Func w_s = Func("w_s");
    w_s(x,y) = 1/ig_vertical(x,y);
    Func w_n = Func("w_n");
    w_n(x,y) = 1/ig_vertical(x,y-2);


    Func w_se = Func("w_se");
    w_se(x,y) = 1/(ig_horizontal(x,y) + ig_vertical(x,y));
    Func w_sw = Func("w_sw");
    w_sw(x,y) = 1/(ig_horizontal(x-2,y) + ig_vertical(x,y));
    Func w_ne = Func("w_ne");
    w_ne(x,y) = 1/(ig_horizontal(x,y) + ig_vertical(x,y-2));
    Func w_nw = Func("w_nw");
    w_nw(x,y) = 1/(ig_horizontal(x-2,y) + ig_vertical(x,y-2));


    Func d_line_top = Func("d_line_top");
    d_line_top(x,y) = green_unenhanced(x,y) - input(x,y);
    Func d_squiggle_top = Func("d_squiggle_top");
    d_squiggle_top(x,y) = (
                                  w_e(x,y) * d_line_top(x+2,y) +
                                  w_w(x,y) * d_line_top(x-2,y) +
                                  w_s(x,y) * d_line_top(x,y+2) +
                                  w_n(x,y) * d_line_top(x,y-2))/(
                                  w_e(x,y) + w_w(x,y) + w_s(x,y) + w_n(x,y)
                          );

    Func d_hat_top_known = Func("d_hat_top_known");
    Expr beta = Expr(0.33f);
    d_hat_top_known(x,y) = beta * d_line_top(x,y) + (1.0f-beta) * d_squiggle_top(x,y);

    Func d_hat_top_diagonal = Func("d_hat_top_diagonal");
    d_hat_top_diagonal(x,y) = (
                                      w_nw(x,y) * d_hat_top_known(x-1,y-1) +
                                      w_ne(x,y) * d_hat_top_known(x+1,y-1) +
                                      w_se(x,y) * d_hat_top_known(x+1,y+1) +
                                      w_sw(x,y) * d_hat_top_known(x-1,y+1))/(
                                      w_nw(x,y) + w_ne(x,y) + w_se(x,y) + w_sw(x,y)
                              );

    Func d_hat_top_g = Func("d_hat_top_g");
    d_hat_top_g(x,y,c) = select(
            x%2 == c, (
                    w_e(x,y) * d_hat_top_diagonal(x+1,y) +
                    w_w(x,y) * d_hat_top_diagonal(x-1,y) +
                    w_s(x,y) * d_hat_top_known(x,y+1) +
                    w_n(x,y) * d_hat_top_known(x,y-1)
            ),
            (
                    w_e(x,y) * d_hat_top_known(x+1,y) +
                    w_w(x,y) * d_hat_top_known(x-1,y) +
                    w_s(x,y) * d_hat_top_diagonal(x,y+1) +
                    w_n(x,y) * d_hat_top_diagonal(x,y-1)
            )
    )/(
                                 w_e(x,y) + w_w(x,y) + w_s(x,y) + w_n(x,y)
                         );

    Func green_debayer = Func("green_debayer");
    green_debayer(x,y) = select(
            x%2+y%2 == 1, input(x,y),
            d_hat_top_known(x,y) + input(x,y));

    //cfapattern = meta['Exif.SubImage1.CFAPattern'].value
    //#import ipdb; ipdb.set_trace()

    Expr redpos, bluepos, red_even, blue_even;
    redpos = 0;
    bluepos = 2;
    red_even = 0;
    blue_even = 1;

    /*if int(cfapattern(0)) == 0:
        redpos = 0
        bluepos = 2
        red_even = 0
        blue_even = 1
    else:
        redpos = 2
        bluepos = 0
        blue_even = True*/

    Func red_debayer = Func("red_debayer");
    red_debayer(x,y) = select(
            x%2+y%2 == redpos, input(x,y),
            x%2+y%2 == bluepos, green_debayer(x,y) - d_hat_top_diagonal(x,y),
            green_debayer(x,y) - d_hat_top_g(x,y, red_even)
    );

    Func blue_debayer = Func("blue_debayer");
    blue_debayer(x,y) = select(
            x%2+y%2 == bluepos, input(x,y),
            x%2+y%2 == redpos, green_debayer(x,y) - d_hat_top_diagonal(x,y),
            green_debayer(x,y) - d_hat_top_g(x,y, blue_even)
    );

    Func debayer = Func("debayer");
    debayer(x,y,c) = select(
            c == 1, green_debayer(x,y),
            c == 0, red_debayer(x,y),
            blue_debayer(x,y)
    );

    Func output = Func("output");
    output(x,y,c) = cast(out_type, debayer(x,y,c));

    output.bound(x, 0, (param.width()/2) * 2);
    output.bound(y, 0, (param.height()/2) * 2);
    output.bound(c, 0, 3);

    Var x_outer, y_outer, x_inner, y_inner, tile_index;
    x_outer = Var("x_outer");
    y_outer = Var("y_outer");
    x_inner = Var("x_inner");
    y_inner = Var("y_inner");
    tile_index = Var("tile_index");

    output.reorder_storage(c, x, y).reorder(c,x,y);
    output.unroll(c);
    //#output.unroll(x,2).unroll(y,2);
    input.compute_root();
    input.unroll(x, 2).unroll(y, 2);
    input.tile(x, y, x_outer, y_outer, x_inner, y_inner, 64, 64);
    input.fuse(x_outer, y_outer, tile_index);
    input.parallel(tile_index);
    input.vectorize(x_inner, 16);

    ig_horizontal.compute_root();
    ig_horizontal.tile(x, y, x_outer, y_outer, x_inner, y_inner, 64, 64);
    ig_horizontal.fuse(x_outer, y_outer, tile_index);
    ig_horizontal.parallel(tile_index);
    ig_horizontal.vectorize(x_inner, 4);

    ig_vertical.compute_root();
    ig_vertical.tile(x, y, x_outer, y_outer, x_inner, y_inner, 64, 64);
    ig_vertical.fuse(x_outer, y_outer, tile_index);
    ig_vertical.parallel(tile_index);
    ig_vertical.vectorize(x_inner, 4);

    pass2.compute_root();
    pass2.tile(x, y, x_outer, y_outer, x_inner, y_inner, 64, 64);
    pass2.fuse(x_outer, y_outer, tile_index);
    pass2.parallel(tile_index);
    pass2.vectorize(x_inner, 4);

    output.compute_root();
    output.tile(x, y, x_outer, y_outer, x_inner, y_inner, 64, 64);
    output.fuse(x_outer, y_outer, tile_index);
    output.parallel(tile_index);
    output.vectorize(x_inner, 4);

    //ig_horizontal.compute_at(output, Var::outermost());
    //ig_vertical.compute_at(output, Var::outermost());

    //green_horizontal.compute_at(pass2, Var::outermost());
    //green_vertical.compute_at(pass2, Var::outermost());

    pmn_h.compute_at(pass2, y_inner); //#.store_at(pass2, y);
    pmn_v.compute_at(pass2, tile_index);
    pmn_d.compute_at(pass2, tile_index);

    //pass2.compute_at(d_hat_top_known, Var::outermost());

    d_hat_top_known.compute_at(output, tile_index);
    //d_line_top.compute_at(d_hat_top_known, Var::outermost());

    //green.compute_at(output, x_inner);
    //output.print_loop_nest();

    //Target target = get_host_target();
    //target.set_feature(Target::Profile);
    //output.compile_jit(target);
    return output;
}