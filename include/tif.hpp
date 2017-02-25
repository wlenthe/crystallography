#ifndef _tif_h_
#define _tif_h_

#include <fstream>
#include <cstdint>
#include <climits>
#include <array>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <vector>

namespace tif {
	//helper struct to hold tif image file directory entries
	struct IfdEntry : public std::array<char, 12> {
		IfdEntry(const std::uint16_t tag = 0x0000, const std::uint16_t type = 0x0000) {
			reinterpret_cast<std::uint16_t*>(this->data())[0] = tag;
			reinterpret_cast<std::uint16_t*>(this->data())[1] = type;
			reinterpret_cast<std::uint32_t*>(this->data())[1] = 0x0001;//count
			reinterpret_cast<std::uint32_t*>(this->data())[2] = 0x0000;//value
		}

		static IfdEntry Byte(const std::uint32_t tag, const std::uint8_t value) {
			IfdEntry ifd(tag, 0x0001);
			reinterpret_cast<std::uint8_t*>(ifd.data()+8)[0] = value;
			return ifd;
		}
		static IfdEntry Short(const std::uint32_t tag, const std::uint16_t value) {
			IfdEntry ifd(tag, 0x0003);
			reinterpret_cast<std::uint16_t*>(ifd.data()+8)[0] = value;
			return ifd;
		}
		static IfdEntry Long(const std::uint32_t tag, const std::uint32_t value) {
			IfdEntry ifd(tag, 0x0004);
			reinterpret_cast<std::uint32_t*>(ifd.data()+8)[0] = value;
			return ifd;
		}

		static IfdEntry ImageWidth(const std::uint32_t width)               {return  Long(0x0100, width);}
		static IfdEntry ImageHeight(const std::uint32_t height)             {return  Long(0x0101, height);}
		static IfdEntry BitsPerSample(const std::uint16_t bit)              {return Short(0x0102, bit);}
		static IfdEntry Compression()                                       {return Short(0x0103, 0x0001);}//no compression
		static IfdEntry PhotometricInterpretation(const std::uint16_t samp) {return Short(0x0106, (3 == samp || 4 == samp) ? 0x0002 : 0x0001);}//rgb / black is zero
		static IfdEntry StripOffset(const std::uint32_t offset)             {return  Long(0x0111, offset);}//black is zero
		static IfdEntry Orientation()                                       {return Short(0x0112, 0x0004);}//origin is bottom left
		static IfdEntry SamplesPerPixel(const std::uint16_t samp)           {return Short(0x0115, samp);}
		static IfdEntry RowsPerStrip(const std::uint32_t rows)              {return  Long(0x0116, rows);}
		static IfdEntry StripByteCount(const std::uint32_t bytes)           {return  Long(0x0117, bytes);}
		static IfdEntry PlanarConfiguration()                               {return Short(0x011c, 0x0001);}//chunky (rgbrgbrgb not rrrgggbbb)
		template <typename T>
		static IfdEntry SampleFormat() {
			std::uint16_t format = 4;//unknown
			if(std::numeric_limits<T>::is_integer) {
				if(!std::numeric_limits<T>::is_signed) format = 1;
				else format = 2;
			} else if(std::numeric_limits<T>::is_iec559) format = 3;
			return Short(0x0153, format);
		}
	};
}

template <typename T>
void writeTif(T const * const data, std::int32_t height, std::int32_t width, std::string fileName, std::uint16_t sampsPerPix = 1) {

	//check if little/big endian and build magic number
	union {
		std::uint32_t i;
		char c[4];
	} u = {0x01020304};
	const bool bigEndian = u.c[0] == 1;
	char bigMagic[4] = {'M','M',0x00,0x2A};
	char litMagic[4] = {'I','I',0x2A,0x00};
	
	//open file and write header
	std::ofstream of(fileName.c_str(), std::ios::out | std::ios::binary);
	std::uint32_t ifdOffset = sizeof(T) * height * width * sampsPerPix + 8;//header, image, ifd
	of.write(bigEndian ? bigMagic : litMagic, 4);
	of.write((char*)&ifdOffset, 4);

	//write data
	// of.write((char*)data.data(), sizeof(T) * data.size());//most readers seam to ignore the orientation flag so data should be written with image convention
	for(int i = height-1; i >= 0; i--) of.write((char*)&data[width * i * sampsPerPix], sizeof(T) * width * sampsPerPix);

	//build ifd entries
	std::vector<tif::IfdEntry> entries;
	entries.push_back(tif::IfdEntry::ImageWidth(width));
	entries.push_back(tif::IfdEntry::ImageHeight(height));
	entries.push_back(tif::IfdEntry::BitsPerSample(CHAR_BIT * sizeof(T)));
	entries.push_back(tif::IfdEntry::Compression());
	entries.push_back(tif::IfdEntry::PhotometricInterpretation(sampsPerPix));
	entries.push_back(tif::IfdEntry::StripOffset(8));
	// entries.push_back(tif::IfdEntry::Orientation());//most readers seam to this flag
	entries.push_back(tif::IfdEntry::SamplesPerPixel(sampsPerPix));
	entries.push_back(tif::IfdEntry::RowsPerStrip(height));
	entries.push_back(tif::IfdEntry::StripByteCount(width * height * sizeof(T) * sampsPerPix));
	entries.push_back(tif::IfdEntry::PlanarConfiguration());
	entries.push_back(tif::IfdEntry::SampleFormat<T>());

	//write IFD
	std::uint16_t numEntries = (std::uint16_t)entries.size();
	of.write((char*)&numEntries, sizeof(numEntries));
	for(size_t i = 0; i < entries.size(); i++) of.write(entries[i].data(), entries[i].size());

	//write offset of next ifd / 0 for last ifd
	std::uint32_t nextOffset = 0;
	of.write((char*)&nextOffset, sizeof(nextOffset));
}

#endif//_tif_h_