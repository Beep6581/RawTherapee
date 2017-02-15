//
// A Halide implementation of "A Low Complexity Color Demosaicing Algorithm Based on Integrated Gradient" by
// King-Hong Chung and Yuk-Hee Chan:
//
// http://www.eie.polyu.edu.hk/~enyhchan/J-JEI-Low_complexity_color_demosaicing_algorithm_based_on_IG.pdf
//

#include <Halide.h>
#include <HalideRuntime.h>

using namespace Halide;

/*
param: Image parameters (mainly type) to compile for.
layout: see below
out_type: the type to cast to at the end

layouts: RGGB, GRBG, GBRG, BGGR
*/

namespace {

class IgdDemosaic : public Halide::Generator<IgdDemosaic> {
public:
    ImageParam in_{Float(32), 2, "input"};
    Param <uint8_t> layout{"layout"};
    Var x, y, c;

    Func build() {
        x = Var("x");
        y = Var("y");
        c = Var("c");
        Expr alpha = 1.0f;

        Func input = Func("input");
        input(x, y) = in_(
                clamp(x, 0, in_.width() - 1),
                clamp(y, 0, in_.height() - 1)
        );

        //#### Horizontal ####
        Func cd_horizontal = Func("cd_horizontal");
        cd_horizontal(x, y) = input(x, y) - (input(x - 1, y) + input(x + 1, y)) / 2;

        Func cd_gradient_horizontal_raw = Func("cd_gradient_horizontal_raw");
        cd_gradient_horizontal_raw(x, y) = (
                                                   absd(cd_horizontal(x, y), cd_horizontal(x + 1, y)) +
                                                   absd(cd_horizontal(x, y), cd_horizontal(x - 1, y))) / 2.0f;

        Func cd_gradient_horizontal = Func("cd_gradient_horizontal");
        cd_gradient_horizontal(x, y) = (cd_gradient_horizontal_raw(x - 1, y) + cd_gradient_horizontal_raw(x, y) +
                                        cd_gradient_horizontal_raw(x + 1, y)) / 3.0f;

        Func ci_gradient_plane_horizontal = Func("ci_gradient_plane_horizontal");
        ci_gradient_plane_horizontal(x, y) = absd(input(x, y), input(x + 2, y));

        Func cd_gradient_plane_horizontal = Func("cd_gradient_plane_horizontal");
        cd_gradient_plane_horizontal(x, y) = cd_gradient_horizontal(x + 1, y);

        Func ig_horizontal = Func("ig_horizontal");
        ig_horizontal(x, y) = ci_gradient_plane_horizontal(x, y) + alpha * (
                2 * cd_gradient_plane_horizontal(x, y) +
                cd_gradient_plane_horizontal(x, y - 1) +
                cd_gradient_plane_horizontal(x, y + 1));

        //#### Vertical ####

        Func cd_vertical = Func("cd_vertical");
        cd_vertical(x, y) = input(x, y) - (input(x, y - 1) + input(x, y + 1)) / 2;

        Func cd_gradient_vertical_raw = Func("cd_gradient_vertical_raw");
        cd_gradient_vertical_raw(x, y) = (
                                                 absd(cd_vertical(x, y), cd_vertical(x, y + 1)) +
                                                 absd(cd_vertical(x, y), cd_vertical(x, y - 1))) / 2.0f;

        Func cd_gradient_vertical = Func("cd_gradient_vertical");
        cd_gradient_vertical(x, y) = (cd_gradient_vertical_raw(x, y - 1) + cd_gradient_vertical_raw(x, y) +
                                      cd_gradient_vertical_raw(x, y + 1)) / 3.0f;

        Func ci_gradient_plane_vertical = Func("ci_gradient_plane_vertical");
        ci_gradient_plane_vertical(x, y) = absd(input(x, y), input(x, y + 2));

        Func cd_gradient_plane_vertical = Func("cd_gradient_plane_vertical");
        cd_gradient_plane_vertical(x, y) = cd_gradient_vertical(x, y + 1);


        Func ig_vertical = Func("ig_vertical");
        ig_vertical(x, y) = ci_gradient_plane_vertical(x, y) + alpha * (
                2 * cd_gradient_plane_vertical(x, y) +
                cd_gradient_plane_vertical(x - 1, y) +
                cd_gradient_plane_vertical(x + 1, y));

        //#### Interpolation ####
        Func green_horizontal = Func("green_horizontal");
        green_horizontal(x, y) = ((input(x - 1, y) + input(x + 1, y)) / 2.0f) +
                                 (2.0f * input(x, y) - input(x - 2, y) - input(x + 2, y)) / 4.0f;

        Func green_vertical = Func("green_vertical");
        green_vertical(x, y) = ((input(x, y - 1) + input(x, y + 1)) / 2.0f) +
                               (2.0f * input(x, y) - input(x, y - 2) - input(x, y + 2)) / 4.0f;

        Func green_neither = Func("green_neither");
        green_neither(x, y) = (green_horizontal(x, y) + green_vertical(x, y)) / 2.0f;


        Func delta_h = Func("delta_h");
        delta_h(x, y) = ig_horizontal(x, y) + ig_horizontal(x - 2, y);

        Func delta_v = Func("delta_v");
        delta_v(x, y) = ig_vertical(x, y) + ig_vertical(x, y - 2);

        Func max_factor = Func("max_factor");
        max_factor(x, y) = max(delta_h(x, y) / delta_v(x, y), delta_v(x, y) / delta_h(x, y));

        float threshold = 1.7f;
        Expr L = 3.0f;


        Func go_to_pass1 = Func("go_to_pass1");
        go_to_pass1(x, y) = select(
                max_factor(x, y) == 1 or max_factor(x, y) > threshold, 1,
                0
        );


        Func pass1 = Func("pass1");
        pass1(x, y) = select(
                delta_h(x, y) == delta_v(x, y), green_neither(x, y),
                delta_h(x, y) > delta_v(x, y), green_vertical(x, y),
                green_horizontal(x, y));


        Func pmn_h = Func("pmn_h");
        pmn_h(x, y) = select(
                go_to_pass1(x, y) == 1, pass1(x, y) - input(x, y),
                green_horizontal(x, y) - input(x, y)
        );


        Func pmn_v = Func("pmn_v");
        pmn_v(x, y) = select(
                go_to_pass1(x, y) == 1, pass1(x, y) - input(x, y),
                green_vertical(x, y) - input(x, y)
        );


        Func pmn_d = Func("pmn_d");
        pmn_d(x, y) = select(
                go_to_pass1(x, y) == 1, pass1(x, y) - input(x, y),
                green_neither(x, y) - input(x, y)
        );

        RDom tDom = RDom(-L, L + 1);

        Func phi_h = Func("phi_h");
        phi_h(x, y) = sum(absd(pmn_h(x, y), pmn_h(x - 2 * tDom.x, y)));

        Func phi_v = Func("phi_v");
        phi_v(x, y) = sum(absd(pmn_v(x, y), pmn_v(x, y - 2 * tDom.x)));

        Func phi_d = Func("phi_d");
        phi_d(x, y) = (
                              sum(absd(pmn_d(x, y), pmn_d(x - 2 * tDom.x, y))) +
                              sum(absd(pmn_d(x, y), pmn_d(x, y - 2 * tDom.x)))) / 2.0f;

        Func smallest_phi = Func("smallest_phi");
        smallest_phi(x, y) = Halide::min(phi_h(x, y), Halide::min(phi_v(x, y), phi_d(x, y)));

        Func pass2 = Func("pass2");
        pass2(x, y) = select(
                go_to_pass1(x, y) == 1, pass1(x, y),
                smallest_phi(x, y) == phi_h(x, y), green_horizontal(x, y),
                smallest_phi(x, y) == phi_v(x, y), green_vertical(x, y),
                green_neither(x, y)
        );

        // Layout definition
        Func is_green = Func("is_green");
        Expr xggb = (layout == 0) || (layout == 3);
        is_green(x, y) = (
                ((x % 2 + y % 2 == 1) && xggb) ||
                ((x % 2 + y % 2 != 1) && !xggb)
        );

        Func is_red = Func("is_red");
        is_red(x, y) = (
                x % 2 + 2 * (y % 2) == layout
        );

        Func is_blue = Func("is_blue");
        is_blue(x, y) = (
                3 - (x % 2 + 2 * (y % 2)) == layout
        );

        Expr red_even = layout % 2 == 0;
        Expr blue_even = !red_even;


        Func green_unenhanced = Func("green_unenhanced");
        green_unenhanced(x, y) = select(
                is_green(x, y), input(x, y),
                pass2(x, y));

        Func w_e = Func("w_e");
        w_e(x, y) = 1 / ig_horizontal(x, y);
        Func w_w = Func("w_w");
        w_w(x, y) = 1 / ig_horizontal(x - 2, y);
        Func w_s = Func("w_s");
        w_s(x, y) = 1 / ig_vertical(x, y);
        Func w_n = Func("w_n");
        w_n(x, y) = 1 / ig_vertical(x, y - 2);


        Func w_se = Func("w_se");
        w_se(x, y) = 1 / (ig_horizontal(x, y) + ig_vertical(x, y));
        Func w_sw = Func("w_sw");
        w_sw(x, y) = 1 / (ig_horizontal(x - 2, y) + ig_vertical(x, y));
        Func w_ne = Func("w_ne");
        w_ne(x, y) = 1 / (ig_horizontal(x, y) + ig_vertical(x, y - 2));
        Func w_nw = Func("w_nw");
        w_nw(x, y) = 1 / (ig_horizontal(x - 2, y) + ig_vertical(x, y - 2));


        Func d_line_top = Func("d_line_top");
        d_line_top(x, y) = green_unenhanced(x, y) - input(x, y);
        Func d_squiggle_top = Func("d_squiggle_top");
        d_squiggle_top(x, y) = (
                                       w_e(x, y) * d_line_top(x + 2, y) +
                                       w_w(x, y) * d_line_top(x - 2, y) +
                                       w_s(x, y) * d_line_top(x, y + 2) +
                                       w_n(x, y) * d_line_top(x, y - 2)) / (
                                       w_e(x, y) + w_w(x, y) + w_s(x, y) + w_n(x, y)
                               );

        Func d_hat_top_known = Func("d_hat_top_known");
        Expr beta = Expr(0.33f);
        d_hat_top_known(x, y) = beta * d_line_top(x, y) + (1.0f - beta) * d_squiggle_top(x, y);

        Func d_hat_top_diagonal = Func("d_hat_top_diagonal");
        d_hat_top_diagonal(x, y) = (
                                           w_nw(x, y) * d_hat_top_known(x - 1, y - 1) +
                                           w_ne(x, y) * d_hat_top_known(x + 1, y - 1) +
                                           w_se(x, y) * d_hat_top_known(x + 1, y + 1) +
                                           w_sw(x, y) * d_hat_top_known(x - 1, y + 1)) / (
                                           w_nw(x, y) + w_ne(x, y) + w_se(x, y) + w_sw(x, y)
                                   );


        Func d_hat_top_g = Func("d_hat_top_g");
        d_hat_top_g(x, y, c) = select(
                // Determine whether known pixels for this colour are above or to the side
                (x % 2 != c), (
                        w_e(x, y) * d_hat_top_diagonal(x + 1, y) +
                        w_w(x, y) * d_hat_top_diagonal(x - 1, y) +
                        w_s(x, y) * d_hat_top_known(x, y + 1) +
                        w_n(x, y) * d_hat_top_known(x, y - 1)
                ),
                (
                        w_e(x, y) * d_hat_top_known(x + 1, y) +
                        w_w(x, y) * d_hat_top_known(x - 1, y) +
                        w_s(x, y) * d_hat_top_diagonal(x, y + 1) +
                        w_n(x, y) * d_hat_top_diagonal(x, y - 1)
                )
        ) / (
                                       w_e(x, y) + w_w(x, y) + w_s(x, y) + w_n(x, y)
                               );

        Func green_debayer = Func("green_debayer");
        green_debayer(x, y) = select(
                is_green(x, y), input(x, y),
                d_hat_top_known(x, y) + input(x, y));

        Func red_debayer = Func("red_debayer");
        red_debayer(x, y) = select(
                is_red(x, y), input(x, y),
                is_blue(x, y), green_debayer(x, y) - d_hat_top_diagonal(x, y),
                green_debayer(x, y) - d_hat_top_g(x, y, red_even)
        );

        Func blue_debayer = Func("blue_debayer");
        blue_debayer(x, y) = select(
                is_blue(x, y), input(x, y),
                is_red(x, y), green_debayer(x, y) - d_hat_top_diagonal(x, y),
                green_debayer(x, y) - d_hat_top_g(x, y, blue_even)
        );

        Func debayer = Func("debayer");
        debayer(x, y, c) = select(
                c == 1, green_debayer(x, y),
                c == 0, red_debayer(x, y),
                blue_debayer(x, y)
        );

        Func output = Func("output");
        output(x, y, c) = debayer(x, y, c);

        output.bound(x, 0, (in_.width() / 2) * 2);
        output.bound(y, 0, (in_.height() / 2) * 2);
        output.bound(c, 0, 3);

        Target target = get_target();
        if(!target.has_gpu_feature()) {
            Var x_outer, y_outer, x_inner, y_inner, tile_index;
            x_outer = Var("x_outer");
            y_outer = Var("y_outer");
            x_inner = Var("x_inner");
            y_inner = Var("y_inner");
            tile_index = Var("tile_index");

            output.reorder_storage(c, x, y).reorder(c, x, y);
            output.unroll(c);

            output.compute_root();
            output.tile(x, y, x_outer, y_outer, x_inner, y_inner, 128, 128);
            output.fuse(x_outer, y_outer, tile_index);
            output.parallel(tile_index);
            output.vectorize(x_inner, 8);

            input.compute_at(output, tile_index);
            input.unroll(x, 2).unroll(y, 2);
            input.vectorize(x, 16);

            ig_horizontal.compute_at(output, tile_index);
            ig_horizontal.vectorize(x, 8);

            cd_gradient_plane_horizontal.compute_at(ig_horizontal, Var::outermost());
            cd_gradient_plane_horizontal.vectorize(x, 16);

            ig_vertical.compute_at(output, tile_index);
            ig_vertical.vectorize(x, 8);

            cd_gradient_plane_vertical.compute_at(ig_vertical, Var::outermost());
            cd_gradient_plane_vertical.vectorize(x, 16);

            pass2.compute_at(output, tile_index);
            pass2.vectorize(x, 8);

            green_horizontal.compute_at(pass2, Var::outermost());
            green_vertical.compute_at(pass2, Var::outermost());

            pmn_h.compute_at(pass2, Var::outermost()).vectorize(x, 8);
            pmn_v.compute_at(pass2, Var::outermost()).vectorize(x, 8);
            pmn_d.compute_at(pass2, Var::outermost()).vectorize(x, 8);

            green_unenhanced.compute_at(output, tile_index);
            green_unenhanced.vectorize(x, 8);

            d_hat_top_known.compute_at(output, tile_index);
            d_hat_top_known.vectorize(x, 8);

            d_hat_top_diagonal.compute_at(output, tile_index);
            d_hat_top_diagonal.vectorize(x, 8);
        } else {
            output.reorder_storage(c, x, y).reorder(c, x, y);
            output.unroll(c);
            //#output.unroll(x,2).unroll(y,2);

            input.compute_root();
            ig_horizontal.compute_root();
            ig_vertical.compute_root();
            pass2.compute_root();
            output.compute_root();

            //pmn_h.compute_at(pass2, Var::outermost());
            //pmn_v.compute_at(pass2, Var::outermost());
            //pmn_d.compute_at(pass2, Var::outermost());
            d_hat_top_known.compute_root();
            d_hat_top_known.vectorize(x, 4);
            d_hat_top_known.gpu_tile(x, y, 8, 8);

            input.unroll(x, 2).unroll(y, 2);
            input.vectorize(x, 4);
            input.gpu_tile(x, y, 8, 8);

            ig_horizontal.vectorize(x, 4);
            ig_horizontal.gpu_tile(x, y, 8, 8);

            ig_vertical.vectorize(y, 4);
            ig_vertical.gpu_tile(x, y, 8, 8);

            pass2.vectorize(x, 4);
            pass2.gpu_tile(x, y, 8, 8);

            output.vectorize(x, 4);
            output.gpu_tile(x, y, 8, 8);
        }

        return output;
    }
};
}

RegisterGenerator<IgdDemosaic> igd_demosaic{"igd_demosaic"};