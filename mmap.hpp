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
#ifndef _mmap_h_
#define _mmap_h_

#include <string>
#include <cstdint>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	#define USE_WIN_MMAP 1

	//limit windows includes
	#ifndef NOMINMAX
		#define NOMINMAX
	#endif
	#ifndef WIN32_LEAN_AND_MEAN
		#define WIN32_LEAN_AND_MEAN
	#endif

	//define windows version (for GetFileSizeEx)
	#ifndef WINVER
		#define WINVER 0x0501
	#endif
	#include <windows.h>
#elif __APPLE__ || __linux__ || __unix__ || defined(_POSIX_VERSION)
	#define USE_POSIX_MMAP 1

	#include <unistd.h>//close
	#include <sys/mman.h>//mmap
	#include <sys/stat.h>
	#include <fcntl.h>//open flags
	#include <errno.h>
	#include <string.h>//strerror
#else
	static_assert(false, "MemoryMappedFile not yet implemented for this platform");
#endif

class MemoryMappedFile {
//helpers to get system error message as a string
#ifdef USE_WIN_MMAP
	static std::string GetErrorMessage() {
		LPSTR buff = NULL;
		DWORD errorCode = GetLastError();
		if(errorCode == 0) return std::string();
		size_t size = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, NULL, errorCode, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR)&buff, 0, NULL);
		std::string message(buff, size);
		LocalFree(buff);
		return message;
	}
	
#elif USE_POSIX_MMAP
	static std::string GetErrorMessage() {return std::string(strerror(errno));}
#endif

	//disable copying
	MemoryMappedFile(MemoryMappedFile const &) = delete;
	void operator=(MemoryMappedFile const &) = delete;

#if USE_WIN_MMAP
	typedef HANDLE FileHandle;
	FileHandle fileHandle;
#elif USE_POSIX_MMAP
	typedef int FileHandle;
#endif
	void* fileBuffer;
	FileHandle fileMap;
	std::string fileName;
	std::uint64_t fileBytes;

	public:
#if USE_WIN_MMAP
		enum class Hint : DWORD {Normal = FILE_ATTRIBUTE_NORMAL, Sequential = FILE_FLAG_SEQUENTIAL_SCAN, Random = FILE_FLAG_RANDOM_ACCESS};
#elif USE_POSIX_MMAP
		enum class Hint : int   {Normal = MADV_NORMAL,           Sequential = MADV_SEQUENTIAL,           Random = MADV_RANDOM};
#endif

		char* rawPointer() {return (char*)fileBuffer;}
		std::uint64_t fileSize() {return fileBytes;}

#if USE_WIN_MMAP
		MemoryMappedFile(std::string filename, Hint hint = Hint::Normal) : fileName(filename), fileHandle(NULL), fileMap(NULL), fileBuffer(NULL) {
			fileHandle = CreateFile(fileName.c_str(), GENERIC_READ, 0, NULL, OPEN_EXISTING, (DWORD)hint, 0);//open file
			if(NULL != fileHandle) {
				LARGE_INTEGER size;
				if(0 != GetFileSizeEx(fileHandle, &size)) {//get file size
					fileBytes = (std::uint64_t)size.QuadPart;
					fileMap = CreateFileMapping(fileHandle, NULL, PAGE_READONLY, 0, 0, NULL);//get pointer to raw data
					if(NULL != fileMap)
						fileBuffer = MapViewOfFile(fileMap, FILE_MAP_READ, 0, 0, 0);//get pointer to raw data
				}
			}
			if(NULL == fileBuffer) {
				if(NULL != fileMap   )	CloseHandle(fileMap   );
				if(NULL != fileHandle)	CloseHandle(fileHandle);
				throw std::runtime_error(fileName + " couldn't be memory mapped: " + GetErrorMessage());
			}
		}

		~MemoryMappedFile() {
			UnmapViewOfFile(fileBuffer);
			CloseHandle(fileMap);
			CloseHandle(fileHandle);
		}

#elif USE_POSIX_MMAP
		MemoryMappedFile(std::string filename, Hint hint = Hint::Normal) : fileName(filename), fileMap(0), fileBuffer(NULL) {
			struct stat64 fileStat;
			if(fstat64(_file, &fileStat) >= 0) {
				fileBytes = fileStat.st_size;
				fileHandle = open(fileName.c_str(), O_RDONLY, 0);//open file
				if(-1 != fileHandle) {
					fileBuffer = mmap64(NULL, fileBytes, PROT_READ, MAP_PRIVATE, fd, 0);
					if(MAP_FAILED == fileBuffer) fileBuffer = NULL;
					else madvise(fileBuffer, fileBuffer, (int)hint | MADV_WILLNEED);
				}
			}
			if(NULL == fileBuffer) {
				if(-1 == fileHandle) close(fd);
				throw std::runtime_error(fileName + " couldn't be memory mapped: " + GetErrorMessage());
			}
		}

		~MemoryMappedFile() {
			munmap(fileBuffer, fileBytes);
			close(fileHandle);
		}
#endif
};

#endif//_mmap_h_