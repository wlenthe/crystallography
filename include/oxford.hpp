/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * Copyright (c) 2018, William C. Lenthe                                           *
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
#ifndef _HKL_H_
#define _HKL_H_

#include <vector>
#include <string>
#include <cstring>
#include <cctype>
#include <sstream>
#include <cstdint>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "symmetry.hpp"
// #include "mmap.hpp"

namespace ctf {
	struct Phase {
		float latticeConstants[6];//a, b, c, alpha, beta, gamma
		std::string name;
		size_t laueGroup, spaceGroup;
	};

	enum class ColumnType {Phase, X, Y, Bands, Error, Euler1, Euler2, Euler3, MAD, BC, BS};

	//templated strtof/strtod/strtold
	template <typename T> T inline strtofp(const char* str, char** str_end) {return std::strtod (str, str_end);}
	template <> float       inline strtofp(const char* str, char** str_end) {return std::strtof (str, str_end);}
#ifdef __MINGW32__
	template <> double      inline strtofp(const char* str, char** str_end) {return (double)std::strtold(str, str_end);}//for some reason strtod is insanely slow in mingw
#else
	template <> double      inline strtofp(const char* str, char** str_end) {return std::strtod (str, str_end);}
#endif
	template <> long double inline strtofp(const char* str, char** str_end) {return std::strtold(str, str_end);}
}

template <typename T>
class HKL {
	public:
		std::string project, author, jobMode;
		size_t xCells, yCells, zCells;
		T xStep, yStep, zStep;
		T acqE1, acqE2, acqE3;
		std::string eulerAngleLine;
		size_t mag, coverage, device;
		T kv, tiltAngle, tiltAxis;
		std::vector<ctf::Phase> phaseList;

		std::vector<T> x, y, error, eulers, mad;
		std::vector<size_t> phases, bands, bc, bs;

		HKL(std::string fileName);
		static bool CanReadFile(std::string fileName) {return CanReadExtension(GetExtension(fileName));}
		std::vector<Symmetry<T> const *> getSymmetries();//convert phase list from laue groups to symmetry pointers

	private:

		static std::string GetExtension(std::string fileName) {
			size_t pos = fileName.find_last_of(".");
			if(std::string::npos == pos) throw std::runtime_error("file name " + fileName + " has no extension");
			std::string ext = fileName.substr(pos+1);
			std::transform(ext.begin(), ext.end(), ext.begin(), [](char c){return std::tolower(c);});
			return ext;
		}

		static bool CanReadExtension(std::string ext) {
			if     (0 == ext.compare("ctf")) return true ;
			else return false;
		}

		size_t readCtf(std::string fileName);
		std::vector<ctf::ColumnType> readCtfHeader(std::istream& is);
		size_t readCtfData(std::istream& is, const std::vector<ctf::ColumnType>& columns);
		size_t readCtfDataMemMap(std::string fileName, std::streamoff offset, const std::vector<ctf::ColumnType>& columns);
};

template <typename T>
HKL<T>::HKL(std::string fileName) {
	size_t pos = fileName.find_last_of(".");
	if(std::string::npos == pos) {
		throw std::runtime_error("file name " + fileName + " has no extension");
	} else {
		//get file extension
		std::string ext = fileName.substr(pos+1);
		std::transform(ext.begin(), ext.end(), ext.begin(), [](char c){return std::tolower(c);});

		//read file according to extension
		size_t pointsRead = 0;
		if(0 == ext.compare("ctf")) pointsRead = readCtf(fileName);
		else throw std::runtime_error("unknown tsl file extension " + ext);

		//check that enough data was read
		if(pointsRead < phases.size()) {
			std::stringstream ss;
			ss << "file ended after reading " << pointsRead << " of " << phases.size() << " data points";
			throw std::runtime_error(ss.str());
		}

		//convert from degrees to radians and 1->0 phase indexing
		const T factor = M_PI / 180.0;
		std::for_each(eulers.begin(), eulers.end(), [factor](T& e){e *= factor;});
		std::for_each(phases.begin(), phases.end(), [](size_t& i){--i;});

		//correct for euler convention

		T* pEulers = eulers.data();
		const T pi = T(2) * std::acos(T(0));
		const T twoPi = T(2) * std::acos(T(0));
		for(size_t i = 0; i < eulers.size() / 3; i++)
			pEulers[3*i] = std::fmod(pEulers[3*i]+pi, twoPi);


	}
}

template <typename T>
std::vector<Symmetry<T> const *> HKL<T>::getSymmetries() {
	std::vector<Symmetry<T> const *> syms(phaseList.size());
	for(size_t i = 0; i < phaseList.size(); i++) {
		switch(phaseList[i].laueGroup) {
			case  1: syms[i] = &symmetry::Groups<T>::Triclinic    ; break;
			case  2: syms[i] = &symmetry::Groups<T>::Monoclinic   ; break;
			case  3: syms[i] = &symmetry::Groups<T>::Orthorhombic ; break;
			case  4: syms[i] = &symmetry::Groups<T>::TrigonalLow  ; break;
			case  5: syms[i] = &symmetry::Groups<T>::Trigonal     ; break;
			case  6: syms[i] = &symmetry::Groups<T>::TetragonalLow; break;
			case  7: syms[i] = &symmetry::Groups<T>::Tetragonal   ; break;
			case  8: syms[i] = &symmetry::Groups<T>::HexagonalLow ; break;
			case  9: syms[i] = &symmetry::Groups<T>::Hexagonal    ; break;
			case 10: syms[i] = &symmetry::Groups<T>::CubicLow     ; break;
			case 11: syms[i] = &symmetry::Groups<T>::Cubic        ; break;
			default: {
				std::stringstream ss;
				ss << "unknown ctf Laue group " << phaseList[i].laueGroup;
				throw std::runtime_error(ss.str());
			}
		}
	}
	return syms;
}

template <typename T>
size_t HKL<T>::readCtf(std::string fileName) {
	//parse header
	std::ifstream is(fileName.c_str());//open file
	if(!is) throw std::runtime_error("ctf file " + fileName + " doesn't exist");
	std::vector<ctf::ColumnType> columnHeaders = readCtfHeader(is);//read header and get column headers

	//allocate memory
	const size_t totalPoints = xCells * yCells * zCells;
	for(const ctf::ColumnType& col : columnHeaders) {
		switch(col) {
			case ctf::ColumnType::Phase : phases.resize(totalPoints    ); break;
			case ctf::ColumnType::X     : x     .resize(totalPoints    ); break;
			case ctf::ColumnType::Y     : y     .resize(totalPoints    ); break;
			case ctf::ColumnType::Bands : bands .resize(totalPoints    ); break;
			case ctf::ColumnType::Error : error .resize(totalPoints    ); break;
			case ctf::ColumnType::Euler1: eulers.resize(totalPoints * 3); break;
			case ctf::ColumnType::Euler2: eulers.resize(totalPoints * 3); break;
			case ctf::ColumnType::Euler3: eulers.resize(totalPoints * 3); break;
			case ctf::ColumnType::MAD   : mad   .resize(totalPoints    ); break;
			case ctf::ColumnType::BC    : bc    .resize(totalPoints    ); break;
			case ctf::ColumnType::BS    : bs    .resize(totalPoints    ); break;
		}
	}

	//read data
	std::streamoff offset = is.tellg();
	is.close();
	try {//attempt to read file with memory map 
		return readCtfDataMemMap(fileName, offset, columnHeaders);
	} catch (...) {//fall back to istream based read
		is.open(fileName.c_str());
		if(!is) throw std::runtime_error("failed to reopen " + fileName + " to read data");
		is.seekg(offset);
		return readCtfData(is, columnHeaders);
	}
}

template <typename T>
std::vector<ctf::ColumnType> HKL<T>::readCtfHeader(std::istream& is) {
	char line[512];
	std::string token;
	bool readProject   = false, readAuthor    = false, readJobMode  = false;
	bool readXCells    = false, readYCells    = false, readZCells   = false;
	bool readXStep     = false, readYStep     = false, readZStep    = false;
	bool readacqE1     = false, readacqE2     = false, readacqE3    = false;
	bool readEulerLine = false;
	bool readMag       = false, readCoverage  = false, readDevice   = false;
	bool readKv        = false, readTiltAngle = false, readTiltAxis = false;
	std::vector<ctf::ColumnType> columnHeaders;

	//check for magic string
	is.getline(line, sizeof(line));//extract first line from file
	if(0 != std::string(line).compare(0, 17, "Channel Text File")) throw std::runtime_error("file is not a valid ctf file (bad header)");//need to compare this way to handle \r\n from windows files on a mac

	//parse header
	while(columnHeaders.empty()) {
		is.getline(line, sizeof(line));//extract entire line from file
		std::istringstream iss(line, sizeof(line));
		if(iss >> token) {//get the key word if the line isn't blank
			//get value for appropriate key
			if       (0 == token.compare("Prj"      )) { project = iss.str(); readProject   = true;
			} else if(0 == token.compare("Author"   )) { author  = iss.str(); readAuthor    = true;
			} else if(0 == token.compare("JobMode"  )) { jobMode = iss.str(); readJobMode   = true;
			} else if(0 == token.compare("XCells"   )) { iss >> xCells      ; readXCells    = true;
			} else if(0 == token.compare("YCells"   )) { iss >> yCells      ; readYCells    = true;
			} else if(0 == token.compare("ZCells"   )) { iss >> zCells      ; readZCells    = true;
			} else if(0 == token.compare("XStep"    )) { iss >> xStep       ; readXStep     = true;
			} else if(0 == token.compare("YStep"    )) { iss >> yStep       ; readYStep     = true;
			} else if(0 == token.compare("ZStep"    )) { iss >> zStep       ; readZStep     = true;
			} else if(0 == token.compare("AcqE1"    )) { iss >> acqE1       ; readacqE1     = true;
			} else if(0 == token.compare("AcqE2"    )) { iss >> acqE1       ; readacqE2     = true;
			} else if(0 == token.compare("AcqE3"    )) { iss >> acqE2       ; readacqE3     = true;
			} else if(0 == token.compare("Euler"    )) {
				eulerAngleLine = std::string(line);
				readEulerLine = true;
				while(iss >> token) {
					if       (0 == token.compare("Mag"      )) { iss >> mag         ; readMag       = true;
					} else if(0 == token.compare("Coverage" )) { iss >> coverage    ; readCoverage  = true;
					} else if(0 == token.compare("Device"   )) { iss >> device      ; readDevice    = true;
					} else if(0 == token.compare("KV"       )) { iss >> kv          ; readKv        = true;
					} else if(0 == token.compare("TiltAngle")) { iss >> tiltAngle   ; readTiltAngle = true;
					} else if(0 == token.compare("TiltAxis" )) { iss >> tiltAxis    ; readTiltAxis  = true;
					}
				}
			} else if(0 == token.compare("Phases"   )) {
				size_t numPhases;
				iss >> numPhases;
				for(size_t i = 0; i < numPhases; i++) {
					//phase lines are formatted as 'a;b;c\talpha;beta;gamma\tname\tlaue_group\tspacegroup'
					is.getline(line, sizeof(line));//extract entire line from file
					std::replace(line, line + sizeof(line), ';', '\t');//convert to fully tab delimited to make things easier
					iss = std::istringstream(line, sizeof(line));

					//sanity check (phase lines should start with a digit)
					iss >> std::skipws;
					if(!std::isdigit(iss.peek())) throw std::runtime_error("end of ctf phases reached before epxected number of phases was read");

					//parse phase
					ctf::Phase p;
					iss >> p.latticeConstants[0] >>  p.latticeConstants[1] >>  p.latticeConstants[2];//a, b, c
					iss >> p.latticeConstants[3] >>  p.latticeConstants[4] >>  p.latticeConstants[5];//alpha, beta, gamma
					iss >> p.name >> p.laueGroup >> p.spaceGroup;
					phaseList.push_back(std::move(p));
				}
			} else if(0 == token.compare("Phase" ) ||
			          0 == token.compare("X"     ) || 
			          0 == token.compare("Y"     ) || 
			          0 == token.compare("Bands" ) || 
			          0 == token.compare("Error" ) || 
			          0 == token.compare("Euler1") || 
			          0 == token.compare("Euler2") || 
			          0 == token.compare("Euler3") || 
			          0 == token.compare("MAD"   ) || 
			          0 == token.compare("BC"    ) || 
			          0 == token.compare("BS"    ) ) {
				bool eulersPresent[3] = {false, false, false};
				do {
					if     (0 == token.compare("Phase" )) columnHeaders.push_back(ctf::ColumnType::Phase );
					else if(0 == token.compare("X"     )) columnHeaders.push_back(ctf::ColumnType::X     );
					else if(0 == token.compare("Y"     )) columnHeaders.push_back(ctf::ColumnType::Y     );
					else if(0 == token.compare("Bands" )) columnHeaders.push_back(ctf::ColumnType::Bands );
					else if(0 == token.compare("Error" )) columnHeaders.push_back(ctf::ColumnType::Error );
					else if(0 == token.compare("Euler1")){columnHeaders.push_back(ctf::ColumnType::Euler1); eulersPresent[0] = true;}
					else if(0 == token.compare("Euler2")){columnHeaders.push_back(ctf::ColumnType::Euler2); eulersPresent[1] = true;}
					else if(0 == token.compare("Euler3")){columnHeaders.push_back(ctf::ColumnType::Euler3); eulersPresent[2] = true;}
					else if(0 == token.compare("MAD"   )) columnHeaders.push_back(ctf::ColumnType::MAD   );
					else if(0 == token.compare("BC"    )) columnHeaders.push_back(ctf::ColumnType::BC    );
					else if(0 == token.compare("BS"    )) columnHeaders.push_back(ctf::ColumnType::BS    );
					else throw std::runtime_error("unknown ctf column header: " + token);
				} while(iss >> token);
				bool allEulers = eulersPresent[0] && eulersPresent[1] && eulersPresent[2];
				bool anyEulers = eulersPresent[0] || eulersPresent[1] || eulersPresent[2];
				if(anyEulers && !allEulers) throw std::runtime_error("only some euler angles present in column headers");
			} else {throw std::runtime_error("unknown ctf header keyword '" + token + "'");}
		}
	}

	//make sure all required parameters were read
	if(!readXCells || !readYCells) throw std::runtime_error("missing ctf dimensions");
	if(!readXStep  || !readYStep ) throw std::runtime_error("missing ctf resolution");
	if(!readacqE1 || !readacqE2 || !readacqE3 ) throw std::runtime_error("missing ctf header AcqE values");
	if(!readEulerLine) throw std::runtime_error("missing ctf euler angle convention line");
	if(!readMag      ) throw std::runtime_error("missing ctf magnification"              );
	if(!readCoverage ) throw std::runtime_error("missing ctf coverage"                   );
	if(!readDevice   ) throw std::runtime_error("missing ctf device"                     );
	if(!readKv       ) throw std::runtime_error("missing ctf accelerating voltage"       );
	if(!readTiltAngle) throw std::runtime_error("missing ctf tilt angle"                 );
	if(!readTiltAxis ) throw std::runtime_error("missing ctf tilt axis"                  );

	//handle missing optional values and return column headers
	if(!readZCells) zCells = 1;
	return columnHeaders;
}

template <typename T>
size_t HKL<T>::readCtfData(std::istream& is, const std::vector<ctf::ColumnType>& columns) {
	//check if the columns are in the normal layout
	const ctf::ColumnType defaultLayout[11] = {
		ctf::ColumnType::Phase ,
		ctf::ColumnType::X     ,
		ctf::ColumnType::Y     ,
		ctf::ColumnType::Bands ,
		ctf::ColumnType::Error ,
		ctf::ColumnType::Euler1,
		ctf::ColumnType::Euler2,
		ctf::ColumnType::Euler3,
		ctf::ColumnType::MAD   ,
		ctf::ColumnType::BC    ,
		ctf::ColumnType::BS    ,
	};
	const bool isDefault = std::equal(columns.cbegin(), columns.cend(), defaultLayout);

	char line[512];
	const size_t totalPoints = xCells * yCells * zCells;

	if(isDefault) {
		for(size_t i = 0; i < totalPoints; i++) {
			char* data = line;
			if(!is.getline(line, sizeof(line))) return i;//get next line
			phases[i    ] = std::strtoul   (data, &data, 10);
			x     [i    ] = ctf::strtofp<T>(data, &data    );
			y     [i    ] = ctf::strtofp<T>(data, &data    );
			bands [i    ] = std::strtoul   (data, &data, 10);
			error [i    ] = ctf::strtofp<T>(data, &data    );
			eulers[3*i  ] = ctf::strtofp<T>(data, &data    );
			eulers[3*i+1] = ctf::strtofp<T>(data, &data    );
			eulers[3*i+2] = ctf::strtofp<T>(data, &data    );
			mad   [i    ] = ctf::strtofp<T>(data, &data    );
			bc    [i    ] = ctf::strtofp<T>(data, &data    );
			bs    [i    ] = ctf::strtofp<T>(data, &data    );
		}
	} else {
		//get a void pointer to each item in order and compute strides
		std::vector<size_t> strides;//size of data type in each column
		std::vector<char* > pointers;//pointers to start of data in each column
		std::vector<size_t> isInt;//1/0 if each column is uint/T
		for(const ctf::ColumnType col : columns) {
			switch(col) {
				case ctf::ColumnType::Phase : strides.push_back(sizeof(phases.front())); pointers.push_back((char*)phases.data()                             ); isInt.push_back(1); break;
				case ctf::ColumnType::X     : strides.push_back(sizeof(x     .front())); pointers.push_back((char*)x     .data()                             ); isInt.push_back(0); break;
				case ctf::ColumnType::Y     : strides.push_back(sizeof(y     .front())); pointers.push_back((char*)y     .data()                             ); isInt.push_back(0); break;
				case ctf::ColumnType::Bands : strides.push_back(sizeof(bands .front())); pointers.push_back((char*)bands .data()                             ); isInt.push_back(1); break;
				case ctf::ColumnType::Error : strides.push_back(sizeof(error .front())); pointers.push_back((char*)error .data()                             ); isInt.push_back(0); break;
				case ctf::ColumnType::Euler1: strides.push_back(sizeof(eulers.front())); pointers.push_back((char*)eulers.data()                             ); isInt.push_back(0); break;
				case ctf::ColumnType::Euler2: strides.push_back(sizeof(eulers.front())); pointers.push_back((char*)eulers.data() + sizeof(eulers.front())    ); isInt.push_back(0); break;
				case ctf::ColumnType::Euler3: strides.push_back(sizeof(eulers.front())); pointers.push_back((char*)eulers.data() + sizeof(eulers.front()) * 2); isInt.push_back(0); break;
				case ctf::ColumnType::MAD   : strides.push_back(sizeof(mad   .front())); pointers.push_back((char*)mad   .data()                             ); isInt.push_back(0); break;
				case ctf::ColumnType::BC    : strides.push_back(sizeof(bc    .front())); pointers.push_back((char*)bc    .data()                             ); isInt.push_back(0); break;
				case ctf::ColumnType::BS    : strides.push_back(sizeof(bs    .front())); pointers.push_back((char*)bs    .data()                             ); isInt.push_back(0); break;
			}
		}

		//parse data
		for(size_t i = 0; i < totalPoints; i++) {
			char* data = line;
			if(!is.getline(line, sizeof(line))) return i;//get next line
			for(size_t j = 0; j < strides.size(); j++) {
				char* pData = pointers[j] + strides[j];
				if(isInt[j])
					reinterpret_cast<size_t*>(pData)[i] = std::strtoul   (data, &data, 10);
				else
					reinterpret_cast<T     *>(pData)[i] = ctf::strtofp<T>(data, &data    );
			}
		}
	}
	return totalPoints;
}

template <typename T>
size_t HKL<T>::readCtfDataMemMap(std::string fileName, std::streamoff offset, const std::vector<ctf::ColumnType>& columns) {
	throw std::runtime_error("not yet implemented for unix");
}

#endif//_HKL_H_