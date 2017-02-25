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
#ifndef _tsl_h_
#define _tsl_h_

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

namespace tsl {
	//structs to store tsl phase data
	struct HKLFamily {
		std::int32_t hkl[3];
		//there is some discrepancy between ang and hdf files for the next 3 terms
		std::int32_t useInIndexing;
		std::int32_t diffractionIntensity;
		std::int32_t showBands;
	};

	struct Phase {
		size_t number;
		std::string materialName, formula, info;
		std::uint32_t symmetry;//tsl symmetry number
		float latticeConstants[6];//a, b, c, alpha, beta, gamma
		std::vector<HKLFamily> hklFamilies;
		float elasticConstants[36];
		std::vector<size_t> categories;
	};

	//enumeration of grid types
	enum class GridType : std::uint8_t {Unknown = 0, Square = 1, Hexagonal = 2};
	std::istream& operator>>(std::istream& is, GridType& grid) {
		std::string name;
		is >> name;
		if     (0 == name.compare("SqrGrid")) grid = GridType::Square;
		else if(0 == name.compare("HexGrid")) grid = GridType::Hexagonal;
		else grid = GridType::Unknown;
		return is;
	}
	std::ostream& operator<<(std::ostream& os, const GridType& grid) {
		switch(grid) {
			case GridType::Square   : return os << "SqrGrid";
			case GridType::Hexagonal: return os << "HexGrid";
			default: throw std::runtime_error("unknown grid type");
		}
	}

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
class TSL {
	public:
		float pixPerUm;
		float xStar, yStar, zStar;
		float workingDistance, xStep, yStep;
		size_t nColsOdd, nColsEven, nRows;
		std::string operatorName, sampleId, scanId;
		tsl::GridType gridType;
		std::vector<tsl::Phase> phaseList;
		std::vector<T> eulers, x, y, iq, ci, sem, fit;
		std::vector<size_t> phases;
		bool tslEulerConvention;

		TSL(std::string fileName);
		static bool CanReadFile(std::string fileName) {return CanReadExtension(GetExtension(fileName));}
		std::vector<Symmetry<T> const *> getSymmetries();//convert phase list from tsl nubmers to symmetry pointers

	private:
		static std::string GetExtension(std::string fileName) {
			size_t pos = fileName.find_last_of(".");
			if(std::string::npos == pos) throw std::runtime_error("file name " + fileName + " has no extension");
			std::string ext = fileName.substr(pos+1);
			std::transform(ext.begin(), ext.end(), ext.begin(), [](char c){return std::tolower(c);});
			return ext;
		}

		static bool CanReadExtension(std::string ext) {
			if     (0 == ext.compare("ang")) return true ;
			else if(0 == ext.compare("osc")) return false;//not quite ready
			else return false;
		}

		void allocate(const size_t tokenCount);
		void transformFromTSL();
		void transformToTSL();

		size_t readAng(std::string fileName);
		size_t readAngHeader(std::istream& is);
		size_t readAngData(std::istream& is, size_t tokens);
		size_t readAngDataMemMap(std::string fileName, std::streamoff offset, size_t tokens);
};

template <typename T>
TSL<T>::TSL(std::string fileName) {
	gridType = tsl::GridType::Unknown;
	size_t pos = fileName.find_last_of(".");
	if(std::string::npos == pos) {
		throw std::runtime_error("file name " + fileName + " has no extension");
	} else {
		//get file extension
		std::string ext = fileName.substr(pos+1);
		std::transform(ext.begin(), ext.end(), ext.begin(), [](char c){return std::tolower(c);});

		//read file according to extension
		size_t pointsRead = 0;
		if(0 == ext.compare("ang")) pointsRead = readAng(fileName);
		else throw std::runtime_error("unknown tsl file extension " + ext);

		//check that enough data was read
		if(pointsRead < iq.size()) {
			std::stringstream ss;
			ss << "file ended after reading " << pointsRead << " of " << iq.size() << " data points";
			throw std::runtime_error(ss.str());
		}

		//correct for tsl orientation reference frame convention
		tslEulerConvention = true;
		transformFromTSL();
	}
}

template <typename T>
std::vector<Symmetry<T> const *> TSL<T>::getSymmetries() {
	std::vector<Symmetry<T> const *> syms(phaseList.size());
	for(size_t i = 0; i < phaseList.size(); i++) {
		switch(phaseList[i].symmetry) {
			case 43: syms[i] = &symmetry::Groups<T>::Cubic        ; break;
			case 23: syms[i] = &symmetry::Groups<T>::CubicLow     ; break;
			case 62: syms[i] = &symmetry::Groups<T>::Hexagonal    ; break;
			case 6 : syms[i] = &symmetry::Groups<T>::HexagonalLow ; break;
			case 42: syms[i] = &symmetry::Groups<T>::Tetragonal   ; break;
			case 4 : syms[i] = &symmetry::Groups<T>::TetragonalLow; break;
			case 32: syms[i] = &symmetry::Groups<T>::Trigonal     ; break;
			case 3 : syms[i] = &symmetry::Groups<T>::TrigonalLow  ; break;
			case 22: syms[i] = &symmetry::Groups<T>::Orthorhombic ; break;
			case 2 : //intentional fall through
			case 20: //intentional fall through
			case 21: syms[i] = &symmetry::Groups<T>::Monoclinic   ; break;
			case 1 : syms[i] = &symmetry::Groups<T>::Triclinic    ; break;
			default: {
				std::stringstream ss;
				ss << "unknown ang symmetry value " << phaseList[i].symmetry;
				throw std::runtime_error(ss.str());
			}
		}
	}
	return syms;
}

template <typename T>
void TSL<T>::allocate(const size_t tokenCount) {
	//compute number of pixels
	size_t totalPoints = 0;
	if(tsl::GridType::Square == gridType) {
		totalPoints = std::max(nColsOdd, nColsEven) * nRows;
	} else if(tsl::GridType::Hexagonal == gridType) {
		totalPoints = size_t(nRows / 2) * (nColsOdd + nColsEven);
		if(1 == nRows % 2) totalPoints += nColsOdd;
	}

	//allocate arrays
	eulers.resize(3*totalPoints);
	x.resize(totalPoints);
	y.resize(totalPoints);
	iq.resize(totalPoints);
	ci.resize(totalPoints);
	phases.resize(totalPoints);
	if(tokenCount > 8) sem.resize(totalPoints);
	if(tokenCount > 9) fit.resize(totalPoints);
}

template <typename T>
void TSL<T>::transformFromTSL() {
	if(tslEulerConvention) {
		tslEulerConvention = false;
		T* pEulers = eulers.data();
		const T pi2 = std::acos(T(0));
		const T twoPi = T(4) * std::acos(T(0));
		for(size_t i = 0; i < eulers.size() / 3; i++)
			pEulers[3*i] = std::fmod(pEulers[3*i]+pi2, twoPi);
	}
}

template <typename T>
void TSL<T>::transformToTSL() {
	if(!tslEulerConvention) {
		tslEulerConvention = true;
		T* pEulers = eulers.data();
		const T pi32  = T(3) * std::acos(T(0));
		const T twoPi = T(4) * std::acos(T(0));
		for(size_t i = 0; i < eulers.size() / 3; i++)
			pEulers[3*i] = std::fmod(pEulers[3*i]+pi32, twoPi);
	}
}

template <typename T>
size_t TSL<T>::readAng(std::string fileName) {
	std::ifstream is(fileName.c_str());//open file
	if(!is) throw std::runtime_error("ang file " + fileName + " doesn't exist");
	size_t tokenCount = readAngHeader(is);//read header and count number of tokens per point
	allocate(tokenCount);

	std::streamoff offset = is.tellg();
	is.close();
	try {//attempt to read file with memory map 
		return readAngDataMemMap(fileName, offset, tokenCount);
	} catch (...) {//fall back to istream based read
		is.open(fileName.c_str());
		if(!is) throw std::runtime_error("failed to reopen " + fileName + " to read data");
		is.seekg(offset);
		return readAngData(is, tokenCount);
	}
}

template <typename T>
size_t TSL<T>::readAngHeader(std::istream& is) {
	char line[512];
	std::string token;
	bool readPixPerUm = false;
	bool readXStar = false, readYStar = false, readZStar = false;
	bool readWorkingDistance = false, readXStep = false, readYStep = false;
	bool readColsOdd = false, readColsEven = false, readRows = false;
	bool readOperatorName = false, readSampleId = false, readScanId = false;
	bool readGridType = false;

	bool readPhaseSymmetry = true;
	size_t targetFamilies = 0, phaseElasticCount = 6;
	bool readPhaseMaterial = true, readPhaseFormula = true, readPhaseInfo = true;
	bool readPhaseLattice = true, readPhaseHkl = true, readPhaseCategories = true;
	while('#' == is.peek()) {//all header lines start with #
		is.getline(line, sizeof(line));//extract entire line from file
		std::istringstream iss(line+1);//skip the '#'
		if(iss >> token) {//get the key word if the line isn't blank
			//get value for appropriate key
			if(0 == token.compare("TEM_PIXperUM")) {iss >> pixPerUm; readPixPerUm = true;
			} else if(0 == token.compare("x-star")) {iss >> xStar; readXStar = true;
			} else if(0 == token.compare("y-star")) {iss >> yStar; readYStar = true;
			} else if(0 == token.compare("z-star")) {iss >> zStar; readZStar = true;
			} else if(0 == token.compare("WorkingDistance")) {iss >> workingDistance; readWorkingDistance = true;
			} else if(0 == token.compare("GRID:")) {iss >> gridType; readGridType = true;
			} else if(0 == token.compare("XSTEP:")) {iss >> xStep; readXStep = true;
			} else if(0 == token.compare("YSTEP:")) {iss >> yStep; readYStep = true;
			} else if(0 == token.compare("NCOLS_ODD:")) {iss >> nColsOdd; readColsOdd = true;
			} else if(0 == token.compare("NCOLS_EVEN:")) {iss >> nColsEven; readColsEven = true;
			} else if(0 == token.compare("NROWS:")) {iss >> nRows; readRows = true;
			} else if(0 == token.compare("OPERATOR:")) {iss >> operatorName; readOperatorName = true;
			} else if(0 == token.compare("SAMPLEID:")) {iss >> sampleId; readSampleId = true;
			} else if(0 == token.compare("SCANID:")) {iss >> scanId; readScanId = true;
			} else if(0 == token.compare("Phase")) {
				//check that all attributes for previous phase were read
				std::stringstream ss;
				ss << phaseList.size();
				if(!readPhaseMaterial) throw std::runtime_error("ang missing material name for phase " + ss.str());
				if(!readPhaseFormula) throw std::runtime_error("ang missing formula for phase " + ss.str());
				if(!readPhaseInfo) throw std::runtime_error("ang missing info for phase " + ss.str());
				if(!readPhaseSymmetry) throw std::runtime_error("ang missing symmetry for phase " + ss.str());
				if(!readPhaseLattice) throw std::runtime_error("ang missing lattice constants for phase " + ss.str());
				if(!readPhaseHkl) throw std::runtime_error("ang missing hkl families for phase " + ss.str());
				if(6 != phaseElasticCount) throw std::runtime_error("ang missing elastic constants for phase " + ss.str());
				if(!readPhaseCategories) throw std::runtime_error("ang missing categories for phase " + ss.str());
				if(!phaseList.empty()) {
					if(targetFamilies < phaseList.back().hklFamilies.size())
						throw std::runtime_error("ang missing some hkl families for phase " + ss.str());
				}
				targetFamilies = phaseElasticCount = 0;
				readPhaseSymmetry = readPhaseMaterial = readPhaseFormula = readPhaseInfo = false;
				readPhaseLattice = readPhaseHkl = readPhaseCategories = false;

				//add a new blank phase to the list
				phaseList.resize(phaseList.size() + 1);
				iss >> phaseList.back().number;
			} else if(0 == token.compare("MaterialName")) {iss >> phaseList.back().materialName; readPhaseMaterial = true;
			} else if(0 == token.compare("Formula")) {iss >> phaseList.back().formula; readPhaseFormula = true;
			} else if(0 == token.compare("Info")) {iss >> phaseList.back().info; readPhaseInfo = true;
			} else if(0 == token.compare("Symmetry")) {iss >> phaseList.back().symmetry; readPhaseSymmetry = true;
			} else if(0 == token.compare("NumberFamilies")) {
				iss >> targetFamilies;
				phaseList.back().hklFamilies.reserve(targetFamilies);
				readPhaseHkl = true;
			} else if(0 == token.compare("LatticeConstants")) {
				iss >> phaseList.back().latticeConstants[0]
				    >> phaseList.back().latticeConstants[1]
				    >> phaseList.back().latticeConstants[2]
				    >> phaseList.back().latticeConstants[3]
				    >> phaseList.back().latticeConstants[4]
				    >> phaseList.back().latticeConstants[5];
				readPhaseLattice = true;
			} else if(0 == token.compare("hklFamilies")) {
				phaseList.back().hklFamilies.resize(phaseList.back().hklFamilies.size() + 1);
				iss >> phaseList.back().hklFamilies.back().hkl[0]
				    >> phaseList.back().hklFamilies.back().hkl[1]
				    >> phaseList.back().hklFamilies.back().hkl[2]
				    >> phaseList.back().hklFamilies.back().useInIndexing
				    >> phaseList.back().hklFamilies.back().diffractionIntensity
				    >> phaseList.back().hklFamilies.back().showBands;
			} else if(0 == token.compare("ElasticConstants")) {
				iss >> phaseList.back().elasticConstants[6*phaseElasticCount + 0]
				    >> phaseList.back().elasticConstants[6*phaseElasticCount + 1]
				    >> phaseList.back().elasticConstants[6*phaseElasticCount + 2]
				    >> phaseList.back().elasticConstants[6*phaseElasticCount + 3]
				    >> phaseList.back().elasticConstants[6*phaseElasticCount + 4]
				    >> phaseList.back().elasticConstants[6*phaseElasticCount + 5];
				++phaseElasticCount;
			} else if(0 == token.compare(0, 10, "Categories")) {//tsl doesn't print space between categories and first number
				//rewind to first character after categories key
				iss.seekg((size_t) iss.tellg() + 10 - token.length());
				std::copy(std::istream_iterator<size_t>(iss), std::istream_iterator<size_t>(), std::back_inserter(phaseList.back().categories));
				readPhaseCategories = true;
			} else {throw std::runtime_error("unknown ang header keyword '" + token + "'");}
		}
	}

	//check that all values of final phase were read
	std::stringstream ss;
	ss << phaseList.size();
	if(!readPhaseMaterial) throw std::runtime_error("ang missing material name for phase " + ss.str());
	if(!readPhaseFormula) throw std::runtime_error("ang missing formula for phase " + ss.str());
	if(!readPhaseInfo) throw std::runtime_error("ang missing info for phase " + ss.str());
	if(!readPhaseSymmetry) throw std::runtime_error("ang missing symmetry for phase " + ss.str());
	if(!readPhaseLattice) throw std::runtime_error("ang missing lattice constants for phase " + ss.str());
	if(!readPhaseHkl) throw std::runtime_error("ang missing hkl families for phase " + ss.str());
	if(6 != phaseElasticCount) throw std::runtime_error("ang missing elastic constants for phase " + ss.str());
	if(!readPhaseCategories) throw std::runtime_error("ang missing categories for phase " + ss.str());
	if(!phaseList.empty()) {
		if(targetFamilies < phaseList.back().hklFamilies.size())
			throw std::runtime_error("ang missing some hkl families for phase " + ss.str());
	}

	//make sure all values were read
	if(!readPixPerUm) throw std::runtime_error("missing ang header value TEM_PIXperUM");
	if(!readXStar) throw std::runtime_error("missing ang header value x-star");
	if(!readYStar) throw std::runtime_error("missing ang header value y-star");
	if(!readZStar) throw std::runtime_error("missing ang header value z-star");
	if(!readWorkingDistance) throw std::runtime_error("missing ang header value WorkingDistance");
	if(!readGridType) throw std::runtime_error("missing ang header value GRID");
	if(!readXStep) throw std::runtime_error("missing ang header value XSTEP");
	if(!readYStep) throw std::runtime_error("missing ang header value YSTEP");
	if(!readColsOdd) throw std::runtime_error("missing ang header value NCOLS_ODD");
	if(!readColsEven) throw std::runtime_error("missing ang header value NCOLS_EVEN");
	if(!readRows) throw std::runtime_error("missing ang header value NROWS");
	if(!readOperatorName) throw std::runtime_error("missing ang header value OPERATOR");
	if(!readSampleId) throw std::runtime_error("missing ang header value SAMPLEID");
	if(!readScanId) throw std::runtime_error("missing ang header value SCANID");

	//extract first line of data without advancing stream
	const std::streamoff start = is.tellg();//save position of data start
	is.getline(line, sizeof(line));//copy first data line
	is.seekg(start);//rewind stream to data start

	//get number of tokens
	size_t tokenCount = 0;
	std::istringstream iss(line);
	while(iss >> token) tokenCount++;
	if(tokenCount < 8) {
		std::stringstream ss;
		ss << "unexpected number of ang values per point (got " << tokenCount << ", expected at least 8)";
		throw std::runtime_error(ss.str());
	}
	return tokenCount;
}

template <typename T>
size_t TSL<T>::readAngData(std::istream& is, size_t tokens) {
	char line[512];//most ang files have 128 byte lines including '\n' so this should be plenty
	bool evenRow = true;
	size_t pointsRead = 0;
	size_t completeRowPoints = 0;
	size_t currentCol = nColsEven - 1;
	const size_t totalPoints = iq.size();
	const bool readSem = tokens > 8;
	const bool readFit = tokens > 9;
	while(pointsRead < totalPoints) {
		char* data = line;
		if(!is.getline(line, sizeof(line))) break;//get next line
		const size_t i = completeRowPoints + currentCol;
		eulers[3*i  ]  = tsl::strtofp<T>(data, &data    );
		eulers[3*i+1]  = tsl::strtofp<T>(data, &data    );
		eulers[3*i+2]  = tsl::strtofp<T>(data, &data    );
		x     [i]      = tsl::strtofp<T>(data, &data    );
		y     [i]      = tsl::strtofp<T>(data, &data    );
		iq    [i]      = tsl::strtofp<T>(data, &data    );
		ci    [i]      = tsl::strtofp<T>(data, &data    );
		phases[i]      = std::strtoul   (data, &data, 10);
		if(readSem) {
			sem[i]     = tsl::strtofp<T>(data, &data    );
			if(readFit) {
				fit[i] = tsl::strtofp<T>(data, NULL     );
			}
		}
		pointsRead++;
		if(0 == currentCol--) {
			completeRowPoints += evenRow ? nColsEven : nColsOdd;
			evenRow = !evenRow;
			currentCol = evenRow ? nColsEven - 1 : nColsOdd - 1;
		}
	}
	return pointsRead;
}

template <typename T>
size_t TSL<T>::readAngDataMemMap(std::string fileName, std::streamoff offset, size_t tokens) {
	throw std::runtime_error("not yet implemented for unix");
/*
	//open memory mapped file
	MemoryMappedFile mapped(fileName, MemoryMappedFile::Hint::Sequential);
	char* data = mapped.rawPointer() + offset;
	std::uint64_t fileBytes = mapped.fileSize();

	//parse data
	bool evenRow = true;
	size_t pointsRead = 0;
	size_t completeRowPoints = 0;
	size_t currentCol = nColsEven - 1;
	size_t totalPoints = iq.size();
	std::uint64_t bytesRead = (std::uint64_t)offset+1;
	const bool readSem = tokens > 8;
	const bool readFit = tokens > 9;
	const bool extraTok = tokens > 10;
	while(pointsRead < totalPoints && bytesRead < fileBytes) {
		char* dStart = data;
		const size_t i = completeRowPoints + currentCol;
		eulers[3*i  ]  = tsl::strtofp<T>(data, &data    );
		eulers[3*i+1]  = tsl::strtofp<T>(data, &data    );
		eulers[3*i+2]  = tsl::strtofp<T>(data, &data    );
		x[i]           = tsl::strtofp<T>(data, &data    );
		y[i]           = tsl::strtofp<T>(data, &data    );
		iq[i]          = tsl::strtofp<T>(data, &data    );
		ci[i]          = tsl::strtofp<T>(data, &data    );
		phases[i]      = std::strtoul   (data, &data, 10);
		if(readSem) {
			sem[i]     = tsl::strtofp<T>(data, &data    );
			if(readFit) {
				fit[i] = tsl::strtofp<T>(data, &data    );
				if(extraTok) while('\n' != *data) ++data;//skip extra tokens until the end of the line
			}
		}
		bytesRead += data - dStart;
		pointsRead++;
		if(0 == currentCol--) {
			completeRowPoints += evenRow ? nColsEven : nColsOdd;
			evenRow = !evenRow;
			currentCol = evenRow ? nColsEven - 1 : nColsOdd - 1;
		}
	}
	return pointsRead;
*/
}

#endif//_tsl_h_