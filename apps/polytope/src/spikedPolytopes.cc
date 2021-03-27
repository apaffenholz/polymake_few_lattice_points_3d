/* Copyright (c) 2020-2021
   Andreas Paffenholz
   Technische Universität Berlin, Germany
   https://www2.mathematik.tu-darmstadt.de/~paffenholz

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
*/

#include "polymake/client.h"
#include "polymake/Array.h"
#include "polymake/Vector.h"
#include "polymake/Matrix.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"

#include "polymake/polytope/normal_form.h"


namespace polymake { namespace polytope {

    Matrix<Integer> homogenize(const Matrix<Integer>& m) {
        Matrix<Integer> ret = ones_vector<Integer>(m.rows())|m;
        return ret;
    }

    std::vector<BigObject> spiked_quasi_minimal_polytopes_impl(int n) {

        std::vector<BigObject> spiked_polytopes;
        char* name = (char*)malloc(100 * sizeof(char));;

        Matrix<Rational> m(4, 4);
        m.col(0) = ones_vector<Rational>(4);

        // 10a
        m[0][1] =  1; m[0][2] =  0;
        m[1][1] =  0; m[1][2] =  2;
        m[2][1] = -1; m[2][2] =  0; m[2][3] = 0;
        m[3][1] =  0; m[3][2] =  0;
        int nred = 2*n-10;
        for ( int a = -1; a <= 0; ++a ) {
            m[0][3] = a;
            for ( int b = -1; b <= 0; ++b ) {
                m[1][3] = b;
                if ( (nred - b) % 3 == 0 && nred - b >= 6 && !( (nred-b)/3 ==2 && a == 0 ) ) {
                    m[3][3] = (nred - b)/3;
                    BigObject p_temp("Polytope<Rational>");
                    p_temp.take("VERTICES") << m;
                    Matrix<Integer> pm = p_temp.call_method("PAIRING_MATRIX");

                    AffineNormalFormTransformation anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm);
                    BigObject p("Polytope<Rational>");
                    p.take("VERTICES") << homogenize(anf.normal_form);
                    sprintf(name,"spiked.p10a.k.%d.a.%d.b.%d",(nred-b)/3,a,b);
                    p.set_name(name);
                    spiked_polytopes.push_back(p);
                }
                if ( (nred - b +1 ) % 3 == 0 && nred - b >= 5 && !( (nred-b+1)/3 ==2 && a == 0 ) ) {
                    m[3][3] = (nred - b + 1)/3;
                    BigObject p_temp("Polytope<Rational>");
                    p_temp.take("VERTICES") << m;
                    Matrix<Integer> pm = p_temp.call_method("PAIRING_MATRIX");

                    AffineNormalFormTransformation anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm);
                    BigObject p("Polytope<Rational>");
                    p.take("VERTICES") << homogenize(anf.normal_form);

                    sprintf(name,"spiked.p10a.k.%d.a.%d.b.%d",(nred-b+1)/3,a,b);
                    p.set_name(name);
                    spiked_polytopes.push_back(p);
                }
            }
        }

        int k = n - 4;
        if (k < 2) {
            return spiked_polytopes;
        }

        // 1
        m[0][1] =  1; m[0][2] = -1; m[0][3] = -1;
        m[1][1] = -1; m[1][2] =  1; m[1][3] =  1;
        m[2][1] = -1; m[2][2] = -1; m[2][3] =  0;
        m[3][1] =  0; m[3][2] =  0; m[3][3] =  k;
        BigObject p1_temp("Polytope<Rational>");
        p1_temp.take("VERTICES") << m;
        Matrix<Integer> pm1 = p1_temp.call_method("PAIRING_MATRIX");

        AffineNormalFormTransformation anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm1);
        BigObject p1("Polytope<Rational>");
        p1.take("VERTICES") << homogenize(anf.normal_form);

        sprintf(name,"spiked.p1.k.%d",k);
        p1.set_name(name);
        spiked_polytopes.push_back(p1);

        // 4
        m[0][1] =  2;
        m[1][2] =  2;
        BigObject p4_temp("Polytope<Rational>");
        p4_temp.take("VERTICES") << m;
        Matrix<Integer> pm4 = p4_temp.call_method("PAIRING_MATRIX");

        anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm4);
        BigObject p4("Polytope<Rational>");
        p4.take("VERTICES") << homogenize(anf.normal_form);

        sprintf(name,"spiked.p4.k.%d",k);
        p4.set_name(name);
        spiked_polytopes.push_back(p4);

        // 5
        m[0][1] =  1;
        m[1][1] =  0; m[1][2] =  1; m[1][3] = -1;
        BigObject p5a_temp("Polytope<Rational>");
        p5a_temp.take("VERTICES") << m;
        Matrix<Integer> pm5a = p5a_temp.call_method("PAIRING_MATRIX");

        anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm5a);
        BigObject p5a("Polytope<Rational>");
        p5a.take("VERTICES") << homogenize(anf.normal_form);

        sprintf(name,"spiked.p5.k.%d.a.-1",k);
        p5a.set_name(name);
        spiked_polytopes.push_back(p5a);

        m[1][3] = 0;
        BigObject p5b_temp("Polytope<Rational>");
        p5b_temp.take("VERTICES") << m;
        Matrix<Integer> pm5b = p5b_temp.call_method("PAIRING_MATRIX");

        anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm5b);
        BigObject p5b("Polytope<Rational>");
        p5b.take("VERTICES") << homogenize(anf.normal_form);
        sprintf(name,"spiked.p5.k.%d.a.",k);
        p5b.set_name(name);
        spiked_polytopes.push_back(p5b);

        // 6
        m[0][2] = 0; m[0][3] = 0;
        m[1][3] = -2;
        BigObject p6a_temp("Polytope<Rational>");
        p6a_temp.take("VERTICES") << m;
        Matrix<Integer> pm6a = p6a_temp.call_method("PAIRING_MATRIX");

        anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm6a);
        BigObject p6a("Polytope<Rational>");
        p6a.take("VERTICES") << homogenize(anf.normal_form);

        sprintf(name,"spiked.p6.k.%d.a.-2",k);
        p6a.set_name(name);
        spiked_polytopes.push_back(p6a);

        m[1][3] = -1;
        BigObject p6b_temp("Polytope<Rational>");
        p6b_temp.take("VERTICES") << m;
        Matrix<Integer> pm6b = p6b_temp.call_method("PAIRING_MATRIX");

        anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm6b);
        BigObject p6b("Polytope<Rational>");
        p6b.take("VERTICES") << homogenize(anf.normal_form);

        sprintf(name,"spiked.p6.k.%d.a.-1",k);
        p6b.set_name(name);
        spiked_polytopes.push_back(p6b);

        m[1][3] = 0;
        BigObject p6c_temp("Polytope<Rational>");
        p6c_temp.take("VERTICES") << m;
        Matrix<Integer> pm6c = p6c_temp.call_method("PAIRING_MATRIX");

        anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm6c);
        BigObject p6c("Polytope<Rational>");
        p6c.take("VERTICES") << homogenize(anf.normal_form);

        sprintf(name,"spiked.p6.k.%d.a.0",k);
        p6c.set_name(name);
        spiked_polytopes.push_back(p6c);

        // 7
        m[0][1] = 2; m[0][2] = 1;
        m[1][1] = -1; m[1][3] = -5;
        BigObject p7a_temp("Polytope<Rational>");
        p7a_temp.take("VERTICES") << m;
        Matrix<Integer> pm7a = p7a_temp.call_method("PAIRING_MATRIX");

        anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm7a);
        BigObject p7a("Polytope<Rational>");
        p7a.take("VERTICES") << homogenize(anf.normal_form);

        sprintf(name,"spiked.p7.k.%d.a.-5",k);
        p7a.set_name(name);
        spiked_polytopes.push_back(p7a);

        m[1][3] = -1;
        BigObject p7b_temp("Polytope<Rational>");
        p7b_temp.take("VERTICES") << m;
        Matrix<Integer> pm7b = p7b_temp.call_method("PAIRING_MATRIX");

        anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm7b);
        BigObject p7b("Polytope<Rational>");
        p7b.take("VERTICES") << homogenize(anf.normal_form);

        sprintf(name,"spiked.p7.k.%d.a.-1",k);
        p7b.set_name(name);
        spiked_polytopes.push_back(p7b);

        k--;
        if ( k < 2 ) {
            return spiked_polytopes;
        }

        // 2
        m[0][1] =  1; m[0][2] = -1; m[0][3] =  0;
        m[1][1] = -1; m[1][2] =  1; m[1][3] = -1;
        m[2][1] = -1; m[2][2] = -1; m[2][3] =  0;
        m[3][1] =  0; m[3][2] =  0; m[3][3] =  k;
        BigObject p2_temp("Polytope<Rational>");
        p2_temp.take("VERTICES") << m;
        Matrix<Integer> pm2 = p2_temp.call_method("PAIRING_MATRIX");

        anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm2);
        BigObject p2("Polytope<Rational>");
        p2.take("VERTICES") << homogenize(anf.normal_form);

        sprintf(name,"spiked.p2.k.%d",k);
        p2.set_name(name);
        spiked_polytopes.push_back(p2);

        // 9
        m[0][1] =  1; m[0][2] =  0; m[0][3] =  0;
        m[1][1] =  0; m[1][2] =  1; m[1][3] =  0;
        m[2][1] = -1; m[2][2] = -1; 
        m[3][1] =  1; m[3][2] =  1;
        for ( int a = -2; a <= 0; ++a ) {
            m[2][3] =  a;
            for ( int b = 0; b <= 1; ++b ) {
                m[3][3] = 2*k-a+b;
                BigObject p9_temp("Polytope<Rational>");
                p9_temp.take("VERTICES") << m;
                Matrix<Integer> pm9 = p9_temp.call_method("PAIRING_MATRIX");

                anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm9);
                BigObject p9("Polytope<Rational>");
                p9.take("VERTICES") << homogenize(anf.normal_form);

                sprintf(name,"spiked.p9.k.%d.a.%d.b.%d",k,a,b);
                p9.set_name(name);
                spiked_polytopes.push_back(p9);                
            }
        }

        // 10b
        m[0][1] =  1; m[0][2] =  0; m[0][3] =  0;
        m[1][1] =  0; m[1][2] =  2;
        m[2][1] = -1; m[2][2] =  0; m[2][2] =  0;
        m[3][1] =  0; m[3][2] =  1; 
        for ( int a = -1; a <= (k == 2 ? -1 : 0); ++a ) {
            m[1][3] = a;
            m[3][3] = k;
            BigObject p10b_temp("Polytope<Rational>");
            p10b_temp.take("VERTICES") << m;
            Matrix<Integer> pm10b = p10b_temp.call_method("PAIRING_MATRIX");

            anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm10b);
            BigObject p10b("Polytope<Rational>");
            p10b.take("VERTICES") << homogenize(anf.normal_form);

            sprintf(name,"spiked.p10b.k.%d.a.%d",k,a);
            p10b.set_name(name);
            spiked_polytopes.push_back(p10b);                
        }

        // 3
        if ( k >= 4 ) {
            m[0][1] =  1; m[0][2] = -1; m[0][3] =  0;
            m[1][1] = -1; m[1][2] =  1; m[0][3] =  0;
            m[2][1] = -1; m[2][2] = -1; m[2][3] =  0;
            m[3][1] =  0; m[3][2] =  0; m[3][3] =  k-1;
            BigObject p3_temp("Polytope<Rational>");
            p3_temp.take("VERTICES") << m;
            Matrix<Integer> pm3 = p3_temp.call_method("PAIRING_MATRIX");

            anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm3);
            BigObject p3("Polytope<Rational>");
            p3.take("VERTICES") << homogenize(anf.normal_form);

            sprintf(name,"spiked.p3.k.%d",k-1);
            p3.set_name(name);
            spiked_polytopes.push_back(p3);
        }

        // 8
        m.resize(5, 4);
        m.col(0) = ones_vector<Rational>(5);
        m[0][1] =  1; m[0][2] =  0; m[0][3] = 0;
        m[1][1] =  0; m[1][2] =  1; m[1][3] = 0;
        m[2][1] = -1; m[2][2] =  0; 
        m[3][1] =  0; m[3][2] = -1; 
        m[4][1] =  0; m[4][2] =  0; m[4][3] = k;
        for ( int a = -1; a <=0; ++a ) {
            for ( int b = (k == 2 && a == 0) ? 3 : a; b < 2*k; ++b ) {
                m[2][3] = a; m[3][3] = b;
                BigObject p8_temp("Polytope<Rational>");
                p8_temp.take("VERTICES") << m;
                Matrix<Integer> pm8 = p8_temp.call_method("PAIRING_MATRIX");

                anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm8);
                BigObject p8("Polytope<Rational>");
                p8.take("VERTICES") << homogenize(anf.normal_form);

                sprintf(name,"spiked.p8.k.%d.a.%d.b.%d",k,a,b);
                p8.set_name(name);
                spiked_polytopes.push_back(p8);
            }
        }

        return spiked_polytopes;
    }

    std::vector<BigObject> spiked_minimal_polytopes_impl(int n)
    {

        std::vector<BigObject> spiked_polytopes;
        char* name = (char*)malloc(100 * sizeof(char));;

        int k = 2 * (n - 5);
        if (k < 4)
        {
            return spiked_polytopes;
        }

        Matrix<Rational> m(4, 4);
        m.col(0) = ones_vector<Rational>(4);
        m[0][1] = 1;
        m[1][2] = 1;
        m[2][1] = -1;
        m[3][2] = -1;
        m[3][3] = k;
        BigObject p1_temp("Polytope<Rational>");
        p1_temp.take("VERTICES") << m;
        Matrix<Integer> pm1 = p1_temp.call_method("PAIRING_MATRIX");

        AffineNormalFormTransformation anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm1);
        BigObject p1("Polytope<Rational>");
        p1.take("VERTICES") << homogenize(anf.normal_form);

        sprintf(name,"spiked.pm.a.0.b.0");
        p1.set_name(name);


        m[3][3] += 1;
        BigObject p2_temp("Polytope<Rational>");
        p2_temp.take("VERTICES") << m;
        Matrix<Integer> pm2 = p2_temp.call_method("PAIRING_MATRIX");

        anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm2);
        BigObject p2("Polytope<Rational>");
        p2.take("VERTICES") << homogenize(anf.normal_form);

        sprintf(name,"spiked.pm.a.0.b.1");
        p2.set_name(name);

        m[2][3] -= 1;
        BigObject p3_temp("Polytope<Rational>");
        p3_temp.take("VERTICES") << m;
        Matrix<Integer> pm3 = p3_temp.call_method("PAIRING_MATRIX");

        anf = affine_normal_form_from_vertices_pairing(dehomogenize(convert_to<Integer>(m)),pm3);
        BigObject p3("Polytope<Rational>");
        p3.take("VERTICES") << homogenize(anf.normal_form);

        sprintf(name,"spiked.pm.a.1.b.1");
        p3.set_name(name);

        spiked_polytopes.push_back(p1);
        spiked_polytopes.push_back(p2);
        spiked_polytopes.push_back(p3);

        return spiked_polytopes;
    }
    

    Array<BigObject> spiked_minimal_polytopes(int n) {

        std::vector<BigObject> spiked_minimal_polytopes;
        if ( n >= 7 ) {
            spiked_minimal_polytopes = spiked_minimal_polytopes_impl(n);
        }

        Array<BigObject> ret(spiked_minimal_polytopes.size());
        for ( unsigned int i = 0; i < spiked_minimal_polytopes.size(); ++i ) {
            ret[i] = spiked_minimal_polytopes[i];
        }

        return ret;
    }

    Array<BigObject> spiked_quasi_minimal_polytopes(int n) {

        std::vector<BigObject> spiked_quasi_minimal_polytopes;
        if ( n >= 7 ) {
            spiked_quasi_minimal_polytopes = spiked_quasi_minimal_polytopes_impl(n);
        }

        Array<BigObject> ret(spiked_quasi_minimal_polytopes.size());
        for ( unsigned int i = 0; i < spiked_quasi_minimal_polytopes.size(); ++i ) {
            ret[i] = spiked_quasi_minimal_polytopes[i];
        }

        return ret;
    }

    Array<BigObject> spiked_polytopes(int n) {

        std::vector<BigObject> spiked_polytopes;
        if ( n >= 7 ) {
            std::vector<BigObject> smp = spiked_minimal_polytopes_impl(n);
            spiked_polytopes.insert(spiked_polytopes.end(), smp.begin(), smp.end());

            std::vector<BigObject> sqmp = spiked_quasi_minimal_polytopes_impl(n);
            spiked_polytopes.insert(spiked_polytopes.end(), sqmp.begin(), sqmp.end());
        }

        Array<BigObject> ret(spiked_polytopes.size());
        for ( unsigned int i = 0; i < spiked_polytopes.size(); ++i ) {
            ret[i] = spiked_polytopes[i];
        }

        return ret;
    }

    UserFunction4perl("# @category Geometry"
                          "# Produce all //3//-dimensional spiked lattice polytopes with //n// lattice points."
                          "# "
                          "# @param Int n the number of lattice points"
                          "# @return Array<Polytope<Rational>>"
                          "# @example This produces all spiked polytopes with 7 lattice points and stores it in the variable $s."
                          "# > $s = spiked_polytopes(7);",
                          &spiked_polytopes, "spiked_polytopes($)");

    UserFunction4perl("# @category Geometry"
                          "# Produce all //3//-dimensional minimal spiked lattice polytopes with //n// lattice points."
                          "# "
                          "# @param Int n the number of lattice points"
                          "# @return Array<Polytope<Rational>>"
                          "# @example This produces all minimal spiked polytopes with 7 lattice points and stores it in the variable $s."
                          "# > $s = spiked_minimal_polytopes(7);",
                          &spiked_minimal_polytopes, "spiked_minimal_polytopes($)");

    UserFunction4perl("# @category Geometry"
                          "# Produce all //3//-dimensional quasi-minimal spiked lattice polytopes with //n// lattice points."
                          "# "
                          "# @param Int n the number of lattice points"
                          "# @return Array<Polytope<Rational>>"
                          "# @example This produces all quasi-minimal spiked polytopes with 7 lattice points and stores it in the variable $s."
                          "# > $s = spiked_quasi_minimal_polytopes(7);",
                          &spiked_quasi_minimal_polytopes, "spiked_quasi_minimal_polytopes($)");

}}