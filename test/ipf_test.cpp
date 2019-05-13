/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * Copyright (c) 2017, William C. Lenthe                                           *
 * All rights reserved.                                                            *
 *                                                                                 *
 * Redistribution and use in source and binary forms, with or without              *
 * modification, are permitted provided that the following conditions are met:     *
 *                                                                                 *
 * 1. Redistributions of source code must retain the above copyright notice, this  *
 *    list of conditions and the following disclaimer.                             *
 *                                                                                 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,    *
 *    this list of conditions and the following disclaimer in the documentation    *
 *    and/or other materials provided with the distribution.                       *
 *                                                                                 *
 * 3. Neither the name of the copyright holder nor the names of its                *
 *    contributors may be used to endorse or promote products derived from         *
 *    this software without specific prior written permission.                     *
 *                                                                                 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"     *
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE       *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE  *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE    *
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL      *
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR      *
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER      *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,   *
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE   *
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <iostream>
#include <fstream>

#include "tif.hpp"
#include "orientation_coloring.hpp"

void hemiLegend(std::string fileName, const size_t dim, void (*ipfFunc)(double const * const, double * const), bool upper = true) {
	double rMax = (double) dim / 2.0;
	std::vector<std::uint8_t> buff(dim * dim * 4, 0x00);//rgba
	for(size_t j = 0; j < dim; j++) {
		double y = (double(j) - rMax) / rMax;
		for(size_t i = 0; i < dim; i++) {
			double x = (double(i) - rMax) / rMax;
			double r = std::hypot(x, y);
			if(r <= 1.0) {
				//unproject to unit sphere (lambert azimuthal equal area projection)
				double n[3], rgb[3];
				double k = 1.0 - x*x - y*y;
				double k2 = std::sqrt(1.0 + k);
				n[0] = k2 * x;
				n[1] = k2 * y;
				n[2] = upper ? k : -k;

				//compute ipf color
				std::uint8_t* pix = &buff[4 * ( (j * dim) + i)];
				ipfFunc(n, rgb);
				pix[0] = (std::uint8_t)std::round(255.0 * rgb[0]);
				pix[1] = (std::uint8_t)std::round(255.0 * rgb[1]);
				pix[2] = (std::uint8_t)std::round(255.0 * rgb[2]);
				pix[3] = 0xFF;//make opaque
			}
		}
	}

	//write image
	writeTif(buff.data(), (int32_t)dim, (int32_t)dim, fileName, 4);
}

int main() {
	try {
		size_t width = 512;
		bool upper = true;

		hemiLegend("-1.tif", width, &coloring::hemiIpf<double>, upper);

		hemiLegend("222.tif", width, &coloring::cyclicIpf<double, 2>, upper);
		hemiLegend("-3.tif", width, &coloring::cyclicIpf<double, 3>, upper);
		hemiLegend("422.tif", width, &coloring::cyclicIpf<double, 4>, upper);
		hemiLegend("622.tif", width, &coloring::cyclicIpf<double, 6>, upper);

		hemiLegend("mmm.tif", width, &coloring::dihedralIpf<double, 2>, upper);
		hemiLegend("-3m1.tif", width, &coloring::dihedralIpf<double, 3>, upper);
		hemiLegend("4_mmm.tif", width, &coloring::dihedralIpf<double, 4>, upper);
		hemiLegend("6_mmm.tif", width, &coloring::dihedralIpf<double, 6>, upper);
		
		hemiLegend("m-3.tif", width, &coloring::cubicLowIpf<double>);
		hemiLegend("m-3m.tif", width, &coloring::cubicIpf<double>);
	} catch (std::exception& e) {
		std::cout << e.what();
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}