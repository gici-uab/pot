// --- this file has a column width of 120 characters ---
/**
 * @mainpage Pairwise Orthogonal Transform
 * 
 * An efficient implementation of the Pairwise Orthogonal Transform by Ian Blanes and Joan Serra-Sagrista.
 * It is a single-thread implementation using a consumer/producer model based on the observer pattern, where
 * the transform is applied line-by-line and minimal memory is required. Both BIL and BSQ modes are supported,
 * although BSQ is only provided for compatibility and as it performs a high amount of seeks file loading performance
 * may be affected.
 * 
 * The transform is described in the article:
 * 
 * I. Blanes and J. Serra-Sagrista, "Pairwise Orthogonal Transform for Spectral Image Coding,"
 * IEEE Transactions on Geoscience and Remote Sensing, vol. 49, no. 3, pp. 961-972, March 2011.
 * http://dx.doi.org/10.1109/TGRS.2010.2071880
 * 
 * Abstract:
 * Spectral transforms are widely used for the codification of remote-sensing imagery, the KLT and wavelets being the
 * two most common transforms. The KLT presents a higher coding performance than wavelets; however, it also carries
 * several disadvantages: high computational cost and memory requirements, a difficult implementation, and its lack of
 * scalability. In this paper, we introduce a novel transform based on the KLT that, while obtaining better coding
 * performance than wavelets, does not have the disadvantages of the KLT. Due to its very small amount of
 * side-information, the transform can be applied in a line-based scheme, which specially reduces the transform memory
 * requirements. Extensive experimental results are conducted for AVIRIS and Hyperion images, both for lossy and
 * lossless, and in combination with various hyperspectral coders. Results for effects on RX anomaly detection and
 * $k$-means clustering are also included. Theoretical and experimental evidence suggest that the proposed transform
 * might be good replacement for wavelets as spectral decorrelator in many of the situations where the KLT is not a
 * suitable option.
 * 
 * Coded by: Ian Blanes.
 * 
 * License:
 * This file is distributed under the terms of the GNU General Public License (GPL) version 3,
 * 
 * If you find it useful, please send an email to any of the authors (so that we can include it in our grant reports).
 * Contact the author for other licensing terms.
 *
 * This work was supported in part by the Spanish Government, by the Catalan Government, and by FEDER under grants
 * TIN2009-14426-C02-01, TIN2009-05737-E/TIN, SGR2009-1224, and FPU2008.
 */
/**
 * @file
 * Pairwise Orthogonal Transform implementation.
 */

#ifndef POT_H
#define POT_H

#include<deque>
#include<vector>
#include<stack>
#include<iostream>
#include<cmath>
#include<cassert>

#include<stdint.h>
//#include<cstdint> when C++0x support is ready

// A half-precision floating point implementation is required.
// Note: Industrial Light & Magic has a one and Mike Acton another (be sure to pick the latest version).
// Mike Acton's is available at http://cellperformance-snippets.googlecode.com/files/half.c
#include "half.h"

/*
 * Side information for each two-component rotation can be reduced by limiting the range of the half-precision
 * floating point number. Negligible quality drops occur.
 */
#define ENABLE_EXPONENT_REDUCTION

/*
 * The real mean of each bandline when coded losslessly is between (-1,1).
 * For images with high bit depth, this mean cannot be assumed to be zero, or an approximate 2dB penalty occurs.
 * Enabling this produces a slight increase in the transform cost.
 */
#define HIGH_PRECISION_COVARIANCE

namespace pot {
using namespace std;

/**
 * Casts from type FROM to short int and saturates the result to fall inside a short int range.
 * Floating point values are properly rounded.
 * @param from input value
 * @param to output value
 */
template <typename FROM> void castWithSaturation (const FROM & from, short & to) {
	if (from > 32767) { // Saturate high
		to = +32767;
	} else if (from < -32768) { // Saturate low
		to = -32768;
	} else {
		to = static_cast<short> (nearbyint(from));
	}
}

/**
 * Casts from type FROM to int. Floating point values are properly rounded.
 * Values are not expected to overflow.
 * @param from input value
 * @param to output value
 */
template <typename FROM> void castWithSaturation (const FROM & from, int & to) {
	to = static_cast<int> (nearbyint(from));

	// while this might be valid, it is very suspicious
	assert(to < 0 || (to & (1<< 30)) == 0);
	assert(to >= 0 || (to & (1<< 30)) != 0);
}

/**
 * Casts from type FROM to short int and saturates the result to fall inside a short int range.
 * Floating point values are properly rounded.
 * @param from input value
 * @param to output value
 */
template <typename FROM> void castWithSaturation (const FROM & from, unsigned short & to) {
	if (from > 65535) { // Saturate high
		to = 65535U;
	} else if (from < 0) { // Saturate low
		to = 0;
	} else {
		to = static_cast<unsigned short> (nearbyint(from));
	}
}

/**
 * Casts from type FROM to int. Floating point values are properly rounded.
 * Values are not expected to overflow.
 * @param from input value
 * @param to output value
 */
template <typename FROM> void castWithSaturation (const FROM & from, unsigned int & to) {
	to = static_cast<unsigned int> (nearbyint(from));

	// while this might be valid, it is very suspicious
	assert(to >= 0 || (to & (1<< 31)) != 0);
}

/**
 * Casts from type FROM to float.
 * @param from input value
 * @param to output value
 */
template <typename FROM> void castWithSaturation (const FROM & from, float & to) {
	to = static_cast<float> (from);

	// No nans introduced here
	assert(! isnan(to) || isnan(from));
}

/**
 * Casts from type FROM to double.
 * @param from input value
 * @param to output value
 */
template <typename FROM> void castWithSaturation (const FROM & from, double & to) {
	to = static_cast<double> (from);

	// No nans introduced here
	assert(! isnan(to) || isnan(from));
}

/**
 * A helper class to handle endianness detection and conversion.
 */
class Endianness {
private:
	/**
	 * Determines whether the current machine is big endian or not. Hopefully this could be optimized by the compiler
	 * at some point in the future (pointer analysis).
	 */
	static bool isBigEndian () {
		union { int i; char j[4]; } a;

		a.i = 0x11000022;
		if (a.j[0] == 0x11) {
			return true; // We are big endian
		} else if (a.j[0] == 0x22) {
			return false; // We are little endian
		} else {
			throw("Mixed endian is not supported.");
		}
	}
	
	/**
	 * Swaps endianness.
	 * @param element input element.
	 * @return a copy of @a element with different endianness.
	 */
	template <typename T> static T swap (const T & element) {
		const size_t elementSize = sizeof(T);
		int nSwaps = elementSize / 2;

		assert (elementSize % 2 == 0);

		union { T t; char c[elementSize]; } a;

		a.t = element;

		for (int i = 0; i < nSwaps; i++) {
			char t = a.c[i];
			a.c[i] = a.c[elementSize - i - 1];
			a.c[elementSize - i - 1] = t;
		}

		return a.t;
	}
	
public:
	/**
	 * Swaps the endianness of @a element as defined in @a fromBigEndian to match the host machine endianness.
	 */
	template <typename T> static T toHostEndianness (const T & element, bool fromBigEndian) {
		bool swap;
		
		if(isBigEndian()) {
			swap = ! fromBigEndian;
		} else {
			swap = fromBigEndian;
		}
		
		if (swap) {
			return Endianness::swap(element);
		} else {
			return element;
		}
	}
	
	/**
	 * Swaps the endianness of an array as defined in @a fromBigEndian to match the host machine endianness.
	 */
	template <typename T> static void toHostEndianness (T * buffer, size_t nElements, bool fromBigEndian) {
		bool swap;
				
		if(isBigEndian()) {
			swap = ! fromBigEndian;
		} else {
			swap = fromBigEndian;
		}
		
		if (swap) {
			for (size_t i = 0; i < nElements; i++) {
				buffer[i] = Endianness::swap(buffer[i]);
			}
		}
	}
};

/**
 * Abstract class for input/output seeking. This class provides offsets for the location of "bandlines" (the portion of
 * a line that contains one band) within an image of specified size. This class has state, and offsets are provided
 * relative to the last location requested and the reading of a bandline at the requested location.
 * 
 * This class is designed to operate in a line by line basis, providing locations within it, and then moving to the
 * next one.
 * 
 * @see BIL and BSQ.
 */
template <typename FileFormatType>
class FileOrder {
protected:
	/// Amount already sought from the next petition.
	streamoff _off;
	const size_t _z, _y, _x;
	
	FileOrder (size_t z, size_t y, size_t x) : _off(0), _z(z), _y(y), _x(x) {}
	
public:
	/**
	 * Provides the offset for the bandline @a index
	 * @param index bandline number. In the range [0, z).
	 */
	virtual streamoff moveTo(int index) = 0;
	
	/**
	 * After being called, offsets will be provided for the next line.
	 */
	virtual void nextLine() = 0;

	/**
	 * Virtual destructor.
	 */
	virtual ~FileOrder() {};
};

/**
 * A derived class from FileOrder for BIL image processing. BIL is Band-Interleaved-by-Line or Row-Interleaved.
 * That means that a image is stored one line at a time and band are stored sequentially within that line.
 * @see FileOrder. 
 */
template <typename FileFormatType>
class BIL : public FileOrder<FileFormatType> {
public:
	BIL (size_t z, size_t y /* can be unknown */, size_t x) : FileOrder<FileFormatType>(z,y,x) {}
	
	streamoff moveTo(int index) {
		streamoff destination = index * this->_x * sizeof(FileFormatType);
		streamoff remaining = destination - this->_off;
		
		this->_off += remaining + this->_x * sizeof(FileFormatType);
		
		return remaining;
	}
	
	void nextLine() {
		this->_off -= this->_z * this->_x * sizeof(FileFormatType);
	}
};

/**
 * A derived class from FileOrder for BSQ image processing. BSQ is Band-Sequential or Band-Interleaved.
 * That means that a image is stored one band at a time and lines are stored sequentially within that band. 
 * @see FileOrder.
 */
template <typename FileFormatType>
class BSQ : public FileOrder<FileFormatType> {
public:
	BSQ (size_t z, size_t y, size_t x) : FileOrder<FileFormatType>(z,y,x) {}
	
	streamoff moveTo(int index) {
		streamoff destination = index * this->_x * this->_y * sizeof(FileFormatType);
		streamoff remaining = destination - this->_off;
		
		this->_off += remaining + this->_x * sizeof(FileFormatType);
		
		return remaining;
	}
	
	void nextLine() {
		this->_off -= this->_x * sizeof(FileFormatType);
	}
};


/**
 * A class to manage transform side information. This class stores side information for one image line and keeps
 * track of the order in which is used. Take into account that side information is read sequentially for each line,
 * but within each line read order is FILO.
 */
template <typename MeanType, typename RotationType>
class SideInformation {
private:
	deque<MeanType> means;
	deque<RotationType> rotations;
	
public:
	/**
	 * Reads the side information required to remove the transform on one line.
	 * @param inputStream the stream where the side information is read.
	 * @param z the amount of bands in a line.
	 */
	void readFromStream (istream & inputStream, const size_t z) {			
		for (size_t i = 0; i < z; i++) {
			MeanType x;
			inputStream.read(reinterpret_cast<char *>(&x), sizeof(MeanType));			
			x = Endianness::toHostEndianness(x, false);
			means.push_back(x);
		}
		
		for (size_t i = 0; i < (z - 1); i++) {
			RotationType x;
			inputStream.read(reinterpret_cast<char *>(&x), sizeof(RotationType));
			x = Endianness::toHostEndianness(x, false);
			rotations.push_back(x);
		}
	}

	/**
	 * Writes the side information produced after applying the transform to one line. Side information written is
	 * discarded.
	 */
	void writeToStream (ostream & outputStream) {
		typename deque<MeanType>::const_iterator i;
		typename deque<RotationType>::const_iterator j;
		
		for (i = means.begin(); i != means.end(); i++) {
			MeanType x = Endianness::toHostEndianness(*i, false);
			outputStream.write(reinterpret_cast<const char *>(&x), sizeof(MeanType));	
		}
		
		for (j = rotations.begin(); j != rotations.end(); j++) {
			RotationType x = Endianness::toHostEndianness(*j, false);
			outputStream.write(reinterpret_cast<const char *>(&x), sizeof(RotationType));
		}
		
		means.clear();
		rotations.clear();  
	}	

	void pushMean(MeanType mean) {
		means.push_back(mean);
	}
	
	MeanType popMean() {
		MeanType mean = means.back();
		means.pop_back();
		return mean;
	}
	
	void pushRotation(RotationType rotation) {
		rotations.push_back(rotation);
	}
	
	RotationType popRotation() {
		RotationType rotation = rotations.back();
		rotations.pop_back();
		return rotation;
	}
};

template <typename T2,typename T1> T1 safeAccumulation (const T2 *const a, const size_t size) {
	T1 * accumulator = new T1[size];

	size_t elements = size;

	for (size_t i = 0; i < elements; i++) {
		accumulator[i] = a[i];
	}

	while (elements > 1) {
		// Add each element to the next. Leave unpaired elements for next round.
		for (size_t i = 0; i < elements / 2; i++) {
			accumulator[i] = accumulator[2 * i] + accumulator[2 * i + 1];
		}

		if (elements % 2 == 0) {
			elements /= 2;
		} else {
			accumulator[elements / 2] = accumulator[elements - 1];
			elements = elements / 2 + 1;
		}
	}

	T1 result = accumulator[0];

	delete accumulator;

	return result;
}

/**
 * Abstract class for a consumer/producer model, similar to the Observer patter. A derived object from this class
 * operates on "bandlines" (the portion of a line that contains one band), which can be sent to it or received from it.
 * 
 * This class also tracks variance of components when are sent to it. By convention negative covariances are invalid.
 * A function is provided to compute covariances if needed. 
 */
template <typename FormatType>
class BandLinePipe {
protected:
	const size_t _x;

	/**
	 * Constructor.
	 * @param x width of the bandlines on which this class operates.
	 */
	BandLinePipe(size_t x) : _x(x) { }
	
	/**
	 * Computes the covariance between to bandlines. Bandlines have zero-mean (because it was forced before), so
	 * the covariance computation is simplified.
	 */
	float computeCovariance (FormatType * bandLine1, FormatType * bandLine2) const {
		double * product = new double[_x];

#ifdef HIGH_PRECISION_COVARIANCE
		double mean1 = safeAccumulation<FormatType,double>(bandLine1, _x) / _x;
		double mean2 = safeAccumulation<FormatType,double>(bandLine2, _x) / _x;

		for (size_t i = 0; i < _x; i++) {
			product[i] = (bandLine1[i] - mean1) * (bandLine2[i] - mean2);
		}
#else
		for (size_t i = 0; i < _x; i++) {
			product[i] = bandLine1[i] * bandLine2[i];
		}
#endif

		double result = safeAccumulation<double,double>(product, _x) / _x;

		delete[] product;

		return result;
	}
	
public:
	/**
	 * Send a bandline to this class.
	 * @param index an index of which of the inputs of this object should receive the bandline.
	 * @param variance if known, the variance of the sent bandline, otherwise -1.
	 */
	virtual void send (FormatType * bandLine, int index, float variance = -1) = 0;
	
	/**
	 * Receive a bandline from this class. No variance tracking is needed in this mode.
	 * @see send().
	 */	
	virtual FormatType * receive (int index) = 0;

	/**
	 * Virtual destructor.
	 */
	virtual ~BandLinePipe() {};
};

/**
 * A derived class from @a BandLinePipe that receives its information from an input stream.
 */
template <typename FormatType, typename FileFormatType>
class BandLineReader : public BandLinePipe<FormatType> {
private:
	istream & _input;
	FileOrder<FileFormatType> & _fileOrderSeeker;
	const bool _bigEndian;
	
public:
	/**
	 * Constructor.
	 * @param inputStream input stream.
	 * @param x width of bandlines (in number element).
	 * @param fileOrderSeeker a @a FileOrder object to set the proper seeking mode within the input stream.
	 * @param bigEndian whether the file is bigEndian or not.
	 */
	BandLineReader(istream & inputStream, size_t x, FileOrder<FileFormatType> & fileOrderSeeker, bool bigEndian) 
		: BandLinePipe<FormatType>(x), _input(inputStream), _fileOrderSeeker(fileOrderSeeker),
		_bigEndian(bigEndian) { }
	
	void send (FormatType * bandLine, int index, float variance) {
		throw("Can't send things to a BandLineReader.");
	}
	
	FormatType * receive (int index) {
		FileFormatType * fileBandLine = new FileFormatType[this->_x];
		FormatType * bandLine = new FormatType[this->_x];
		
		_input.seekg(_fileOrderSeeker.moveTo(index), ios_base::cur);
		_input.read(reinterpret_cast<char *>(fileBandLine), this->_x * sizeof(FileFormatType));
		
		Endianness::toHostEndianness<FileFormatType> (fileBandLine, this->_x, _bigEndian);		
		for (size_t i = 0; i < this->_x; i++) {
			castWithSaturation(fileBandLine[i], bandLine[i]);
		}
		
		delete[] fileBandLine;
		
		return bandLine;
	}
	
	/**
	 * Indicate to this object that all bands for a line have been already requested and that it should move its
	 * operations to provide bandlines for the next line.
	 */
	void startNextLine() {
		_fileOrderSeeker.nextLine();
	}
};

/**
 * A derived class from @a BandLinePipe that sends its information to an output stream.
 * @see BandLineReader.
 */
template <typename FormatType, typename FileFormatType>
class BandLineWriter : public BandLinePipe<FormatType> {
	
private:
	ostream & _output;
	FileOrder<FileFormatType> & _fileOrderSeeker;
	const bool _bigEndian;
	
public:
	BandLineWriter(ostream & outputStream, size_t x, FileOrder<FileFormatType> & fileOrderSeeker, bool bigEndian) 
		: BandLinePipe<FormatType>(x), _output(outputStream), _fileOrderSeeker(fileOrderSeeker), 
		_bigEndian(bigEndian) { }
	
	void send (FormatType * bandLine, int index, float variance = -1) {
		FileFormatType * fileBandLine = new FileFormatType[this->_x];
		
		for (size_t i = 0; i < this->_x; i++) {
			castWithSaturation(bandLine[i], fileBandLine[i]);
		}
		
		// toHostEndianness is an involutary function
		Endianness::toHostEndianness<FileFormatType> (fileBandLine, this->_x, _bigEndian);
		
		_output.seekp(_fileOrderSeeker.moveTo(index), ios_base::cur);	
		_output.write(reinterpret_cast<char *>(fileBandLine), this->_x * sizeof(FileFormatType));
		
		delete[] bandLine;
		delete[] fileBandLine;
	}
	
	FormatType * receive (int index) {
		throw("Can't receive things from a BandLineReader.");
	}
	
	void startNextLine() {
		_fileOrderSeeker.nextLine();
	}
};

/**
 * An abstract class representing one of the pair transforms. This class is abstract so that its virtual methods can be
 * overridden with proper lossy and lossless variants of the transform.
 * 
 * Pairs of bandlines can be either sent or received to/from this class. When two bandlines are sent to it (index = 0
 * and index = 1), it applies a transform to them, and sends them to the outputs specified in the class constructor.
 * When a bandline is received from it, this class requests a pair of bandlines from its outputs, removes a
 * transform from them (unless the results is caches for the request of the other bandline of a pair), and sends
 * one of them as the requested result.
 */
template <typename FormatType, typename MeanType, typename RotationType>
class PairTransform : public BandLinePipe<FormatType> {
private:
	BandLinePipe<FormatType> & _output1, & _output2;
	const int _index1, _index2;
	
	FormatType * _bandLine1, * _bandLine2;
	float _variance1, _variance2;
	
	/**
	 * Copy constructor. Only the original copy owns any bandline it already had.
	 */
	PairTransform (const PairTransform& a) { _bandLine1 = 0; _bandLine2 = 0; };

	/**
	 * Assignament operator. Only the original copy owns any bandline it already had.
	 */
	PairTransform& operator=(const PairTransform& a) { _bandLine1 = 0; _bandLine2 = 0; };

protected:
	SideInformation<MeanType, RotationType> & _sideInformation;

	/**
	 * Applies a pair transform to @a bandLine1 and @a bandLine2. Variances are already computed and passed as
	 * parameter.
	 * @return the resulting variance of the first bandline (the second is not needed and is discarded). 
	 */
	virtual float apply (FormatType * bandLine1, float variance1, FormatType * bandLine2, float variance2) = 0;
	
	/**
	 * Does the opposite of @a apply().
	 */
	virtual void remove (FormatType * bandLine1, FormatType * bandLine2) = 0;
	
public:
	/**
	 * Constructor.
	 * @param output1 output where the first result of a pair transform is sent.
	 * @param index1 index associated with output1, where results are sent.
	 * @param output2 idem.
	 * @param index2 idem.
	 * @param x bandline width.
	 * @param sideInformation object where sideInformation is accumulated (this object is shared).
	 */
	PairTransform (BandLinePipe<FormatType> & output1, int index1, BandLinePipe<FormatType> & output2, int index2,
		size_t x, SideInformation<MeanType, RotationType> & sideInformation)
		: BandLinePipe<FormatType>(x), _output1(output1), _output2(output2), _index1(index1), _index2(index2),
		_bandLine1(0), _bandLine2(0), _variance1(0), _variance2(0), _sideInformation(sideInformation) { }

	void send (FormatType * bandLine, int index, float variance = -1) {		
		if (variance < 0) 
			variance = this->computeCovariance(bandLine,bandLine);
			
		if (index == 0) {
			_bandLine1 = bandLine; 
			_variance1 = variance;
		} else if (index == 1) {
			_bandLine2 = bandLine;
			_variance2 = variance;
		}
		
		if (_bandLine1 && _bandLine2) {
			float variance = apply (_bandLine1, _variance1, _bandLine2, _variance2);
			
			_output1.send(_bandLine1, _index1, variance);
			_output2.send(_bandLine2, _index2);
			// we are not the owners what we have already sent 
			_bandLine1 = 0;
			_bandLine2 = 0;
		}
	}
	
	FormatType * receive (int index) {
		if (! _bandLine2 && ! _bandLine1) {
			// Once both of them have been requested, compute a new pair
			_bandLine2 = _output2.receive(_index2);
			_bandLine1 = _output1.receive(_index1);
			remove (_bandLine1, _bandLine2);
		}

		if (index == 1) {
			FormatType * t = _bandLine2;
			_bandLine2 = 0;
			return t;
		} else /*if (index == 0)*/ {
			FormatType * t = _bandLine1;
			_bandLine1 = 0;
			return t;
		}
	}

	/**
	 * Virtual destructor.
	 */
	virtual ~PairTransform() {};
};

/**
 * A half-float conversion class that warps any of the available conversion code.
 * In this case Mike Acton's one.
 */
class HalfFloat {
	union FLOAT_32 { float f32; uint32_t u32; };

public:
	static float halfToFloat(uint16_t a) {
		union FLOAT_32 out;
		out.u32 = half_to_float(a);
		return out.f32;
	}

	static uint16_t floatToHalf(float a) {
		union FLOAT_32 in;
		in.f32 = a;
		return half_from_float(in.u32);
	}

};

/**
 * Lossy specialization of @a PairTransform. While all the other code is just how to perform the transform in a very
 * efficient line-based mode, this is the real deal of the transform.
 * @see PairTransform.
 */
template <typename FormatType, typename MeanType, typename RotationType>
class LossyPairTransform : public PairTransform<FormatType, MeanType, RotationType> {	
protected:	
	float apply (FormatType * bandLine1, float a, FormatType * bandLine2, float d) {
		
		// Compute the KLT for bandLine1 and bandLine2
		float t;
		
		float b = this->computeCovariance(bandLine1, bandLine2);
		float s = sqrtf((a-d) * (a-d) + 4.0 * b * b); 
		float l1;
	
		assert(! isnan(s) && ! isnan(b));
	
		const float SMALL_S_THRESHOLD = 0.0001f;
		if (s > SMALL_S_THRESHOLD) {
			t /* = -q */ = copysignf(sqrtf(abs(0.5 - (a-d) / (2.0 * s))), b); // zero has sign here
			// An abs() is here because (0.5 - (a-d) / (2.0 * s)) might be between [-0.001, 0)
		} else {
			t = 0.0;
		}

#ifdef ENABLE_EXPONENT_REDUCTION
		// Reduce exponent pressure (by 7 bits)
		if (abs(t) < 0.00390625f) {
			t = +1.0f * 0.0;
		} else if (abs(t) < 0.0078125f) {
			t = copysign(0.0078125f, t);
		}
#endif

		assert (! isnan(t));

		float oldt = t;

		uint16_t store = HalfFloat::floatToHalf(t);
		t = HalfFloat::halfToFloat(store);

		assert (! isnan(t));

		// Assert a sign stable conversion
		assert(t * oldt >= 0);

		// Make sure we are not sending negative zeros
		if (t == 0) {
			t = +1.0f * 0.0;
		}

		// Assert that the conversion does not remove more precision than it has to
		if (abs(t - oldt) > 0.01) {
			cerr << "error: " << t << " " << oldt << endl;
			assert(0);
		}
		
		float p /* = u */ = sqrtf(1.0 - t * t);
		
		// Apply new transform
		// x' = px + ty, y' = qx + uy;
		
		for (size_t i = 0; i < this->_x; i++) {
			FormatType x = bandLine1[i], y = bandLine2[i];
				
			FormatType xp =  p * x + t * y;
			FormatType yp = -t * x + p * y;
				
			bandLine1[i] = xp;
			bandLine2[i] = yp;
		}
		
		// Compute the variance of bandLine1 (lambda_1 from the KLT).
		l1 = (a + d + s) / 2.0;
		
		// Only one of the variables needs to be stored as side information.
		this->_sideInformation.pushRotation(store);
		return l1;
	}
	
	void remove (FormatType * bandLine1, FormatType * bandLine2) {
		RotationType store = this->_sideInformation.popRotation();
		float t /* = -q */ = HalfFloat::halfToFloat(store);
		
		float p /* = u */ = sqrtf(1.0 - t * t);

		
		// Remove transform
		// x' = px + qy, y' = tx + uy;
		
		for (size_t i = 0; i < this->_x; i++) {
			FormatType x = bandLine1[i], y = bandLine2[i];
				
			FormatType xp = p * x - t * y;
			FormatType yp = t * x + p * y;

			bandLine1[i] = xp;
			bandLine2[i] = yp;
		}
	}
public:
	LossyPairTransform (BandLinePipe<FormatType> & output1, int index1, BandLinePipe<FormatType> & output2, int index2,
		size_t x, SideInformation<MeanType, RotationType> & sideInformation)
		: PairTransform<FormatType, MeanType, RotationType> (output1, index1, output2, index2,
		x, sideInformation) {}
};

/**
 * Lossless specialization of @a PairTransform. Basically is the same as @a LossyPairTransform, but with an extra
 * lifting stage for perfect recovery.
 * @see PairTransform, LossyPairTransform.
 */
template <typename FormatType, typename MeanType, typename RotationType>
class LosslessPairTransform : public PairTransform<FormatType, MeanType, RotationType>  {
private:

	/**
	 * Decompose a transform into lifting steps.
	 * This has to be strictly the same in the application and removal.
	 */
	void lifting (const float t, double & w1, double & w2, bool & permutationRequired) const {
		// Lifting decomposition
		double p /* = u */ = sqrt(1.0 - (double)t * t);
		
		if (abs(t) >= p) {
			w1 = (p - 1.0) / t;
			w2 = t;
			//w3 = (p - 1) / t;
			
			assert (abs(w1) < 1.1); // 1 + margin
			assert (abs(w2) < 1.1); // 1 + margin

			permutationRequired = false;
		} else {
			w1 = (1.0 - t) / p;
			w2 = -p;
			//w3 = (1 - t) / p;
			
			assert (w1 < 2.51 && w1 > 0.40); // 1 + sqrt(2) + margin, sqrt(2) - 1 - margin
			assert (abs(w2) < 1.1); // 1 + margin

			permutationRequired = true;
		}
	}
	
	/**
	 * Adds to FormatType values and checks for possible overflows.
	 * result = a + sign*[w*b]
	 */
	FormatType LiftingStep (FormatType a, int sign, double w, FormatType b) {
		FormatType c = (sign >= 0? 1: -1) * lrint (w * b);

		assert (a <= 0 || a + c >= c);
		assert (a >= 0 || a + c <= c);

		return a + c;
	}


protected:
	float apply (FormatType * bandLine1, float a, FormatType * bandLine2, float d) {
		
		// Compute the KLT for bandLine1 and bandLine2
		float t;
		
		float b = this->computeCovariance(bandLine1, bandLine2);
		float s = sqrtf((a-d) * (a-d) + 4.0 * b * b); 
		float l1;
		
		const float SMALL_S_THRESHOLD = 0.0001f;
		if (s > SMALL_S_THRESHOLD) {
			t /* = -q */ = copysignf(sqrtf(abs(0.5 - (a-d) / (2.0 * s))), b); // zero has sign here
			// An abs() is here because (0.5 - (a-d) / (2.0 * s)) might be between [-0.001, 0)
		} else {
			t = 0.0;
		}

#ifdef ENABLE_EXPONENT_REDUCTION
		// Reduce exponent pressure (by 7 bits)
		if (abs(t) < 0.00390625f) {
			t = +1.0f * 0.0;
		} else if (abs(t) < 0.0078125f) {
			t = copysign(0.0078125f, t);
		}
#endif

		// Convert to a half float for storage
		assert (! isnan(t));

		uint16_t store = HalfFloat::floatToHalf(t);
		t = HalfFloat::halfToFloat(store);
		
		assert (! isnan(t));

		// Make sure we are not sending negative zeros
		if (t == 0) {
			t = +1.0f * 0.0;
		}
		
		// Lifting decomposition
		double w1, w2;
		bool permutationRequired;
		
		lifting (t, w1, w2, permutationRequired);
		
		// Apply new transform
		// y += [w1 * x]; x += [w2 * y]; y += [w3 * x]; 
		
		for (size_t i = 0; i < this->_x; i++) {
			FormatType x = bandLine1[i], y = bandLine2[i];
			
			y = LiftingStep (y, 1, w1, x);
			x = LiftingStep (x, 1, w2, y);
			y = LiftingStep (y, 1, w1, x);
			
			if (! permutationRequired) {
				bandLine1[i] = x;
				bandLine2[i] = y;
			} else {
				bandLine1[i] = y;
				bandLine2[i] = x;
			}
		}
		
		// Compute the variance of bandLine1 (lambda_1 from the KLT).
		l1 = (a + d + s) / 2.0;
		
		// Only one of the variables needs to be stored as side information.
		this->_sideInformation.pushRotation(store);
		return l1;
	}
	
	void remove (FormatType * bandLine1, FormatType * bandLine2) {
		RotationType store = this->_sideInformation.popRotation();
		float t /* = -q */ = HalfFloat::halfToFloat(store);
		
		// Lifting decomposition
		double w1, w2;
		bool permutationRequired;
		
		lifting (t, w1, w2, permutationRequired);
		
		// Remove transform
		// y -= [w3 * x]; x -= [w2 * y]; y -= [w1 * x]; 
		
		for (size_t i = 0; i < this->_x; i++) {
			FormatType x, y;
			
			if (! permutationRequired) {
				x = bandLine1[i];
				y = bandLine2[i];
			} else {
				y = bandLine1[i];
				x = bandLine2[i];
			}

			y = LiftingStep (y, -1, w1, x);
			x = LiftingStep (x, -1, w2, y);
			y = LiftingStep (y, -1, w1, x);
				
			bandLine1[i] = x;
			bandLine2[i] = y;
		}
	}
public:	
	LosslessPairTransform (BandLinePipe<FormatType> & output1, int index1, BandLinePipe<FormatType> & output2, 
		int index2, size_t x, SideInformation<MeanType, RotationType> & sideInformation)
		: PairTransform<FormatType, MeanType, RotationType> (output1, index1, output2, index2, x, sideInformation) {}
};


/**
 * Main class of the transform. This class is the one that builds the composition of pairwise transforms and applies it
 * to images.
 */
template <typename FormatType, typename FileFormatType, typename MeanType, typename RotationType>
class PairwiseOrthogonalTransform {
private:
	size_t _x, _y, _z;
	bool _bigEndian;
	bool _isBSQ;
	bool _isLossless;
	
	vector<PairTransform<FormatType, MeanType, RotationType> *> pairs;
	SideInformation<MeanType, RotationType> sideInformation;
	
	/**
	 * Creates the composition of pair transforms. 
	 * @param terminalElement the element where the second bandline of each pair will go after transformed.
	 */
	void createStructure(BandLinePipe<FormatType> * const terminalElement) {
		
		// Calculate the number of components remaining at each level.
		size_t components = _z;		
		vector<size_t> s;
		
		while (components > 1) {
			s.push_back(components);
			components = (components / 2) + (components % 2);
		}
		
		// first input from terminalElement not used.
		int currentFileIndex = 0;
		
		// Components up to two levels upper are required for levels with an odd number of components.
		int currentLevelStart = 0;
		int upperLevelStart = 0;
		int upperUpperLevelStart = 0;
		
		// For each level starting at the last one		
		for (int j = s.size() - 1; j >= 0; j--) { 
			size_t componentsThisLevel = s[j];
			
			// Add clusters as needed
			for (size_t i = 0; i < componentsThisLevel / 2; i++) {
				BandLinePipe<FormatType> * output1, * output2;
				int index1, index2;
				
				// Determine this clusters outputs
				if (componentsThisLevel == 2) {
					/* top level: both go to the BandWriter */
					output1 = output2 = terminalElement;
					index1 = currentFileIndex++;
					index2 = currentFileIndex++;
				} else {
					/* remaining levels */
					int outputPosition = i;
										
					if (componentsThisLevel % 2 == 1 && j % 2 == 1) {
						// We have one remaining component in this level
						// and it is aligned at the left side
						outputPosition++;
					}
					
					// Careful here with the unpaired component (on the first level should be aligned right).		
					int outputPairIndex;
					int componentsUpperLevel = s[j+1];
					
					if (componentsUpperLevel % 2 == 1 && j % 2 == 1 && outputPosition == componentsUpperLevel - 1) {						
						// Our output component corresponds with the forwarded component of the next level (right side)
						outputPairIndex = upperLevelStart - 1;
						index1 = 1;
					} else if (componentsUpperLevel % 2 == 1 && j % 2 == 0 && outputPosition == 0) {
						// Our output component corresponds with the forwarded component of the next level (left side)
						outputPairIndex = upperUpperLevelStart;
						index1 = 0;
					} else {
						if (componentsUpperLevel % 2 == 1 && j % 2 == 0) {
							// If in the next level the first component is pushed forward
							// the relative position of this input is decreased by one
							outputPosition--;
						}
						
						outputPairIndex = upperLevelStart + (outputPosition) / 2;
						index1 = (outputPosition) % 2;
					}
					
					output1 = pairs[outputPairIndex]; 
					
					// The second output of a pair transform always goes to the terminalElement.
					output2 = terminalElement;
					index2 = currentFileIndex++;
				}
				
				// Once the outputs are determined, create the pair transform.
				PairTransform<FormatType, MeanType, RotationType> * pairTransform;
				
				if (_isLossless) {
					pairTransform = new LosslessPairTransform<FormatType, MeanType, RotationType>(*output1, index1,
						*output2, index2, _x, sideInformation);
				} else {
					pairTransform = new LossyPairTransform<FormatType, MeanType, RotationType>(*output1, index1,
						*output2, index2, _x, sideInformation);
				}
					
				pairs.push_back(pairTransform);
			}
			
			// Next level
			upperUpperLevelStart = upperLevelStart;
			upperLevelStart = currentLevelStart;
			currentLevelStart += componentsThisLevel / 2;
		}
	}
	
	/**
	 * Destroy a pair structure.
	 */
	void destroyStructure() {
		typename vector<PairTransform<FormatType, MeanType, RotationType> *>::iterator i;
		
		for (i = pairs.begin(); i != pairs.end(); i++) {
			delete *i;
		}
		
		pairs.clear();
	}
	
	/**
	 * Shifts the elements of a bandline so that is has zero mean.
	 */
	void applyZeroMean (FormatType * bandLine) {
		// Accumulation is in a double because otherwise for some
		// data combinations we may have an overflow
		double total = safeAccumulation<FormatType, double>(bandLine, _x);

		// Rounding might be needed here for increased accuracy,
		// but experimental results don't show any SNR difference.
		// MeanType m = rint(total / _x);
		MeanType m;
		castWithSaturation(total / (double) _x, m);
		
		for (size_t i = 0; i < _x; i++) {
			//bandLine[i] = rint(bandLine[i] - m);
			bandLine[i] -= m;
		}
		
		this->sideInformation.pushMean(m);
	}
	
	/**
	 * Reverses the application of @a applyZeroMean().
	 */
	void removeZeroMean (FormatType * bandLine) {
		MeanType m = this->sideInformation.popMean();
		
		for (size_t i = 0; i < _x; i++) {
			//bandLine[i] = rint(bandLine[i] + m);
			bandLine[i] += m;
		}
	}
	
	/**
	 * A class that forces all exceptions on a stream. This class saves the previous stream exception flags
	 * and restores them at destruction, so that only while the object is alive exceptions are forced.  
	 */
	class ForceExceptions {
		ios & _o;
		ios_base::iostate _s;
	public:
		ForceExceptions(ios & o) : _o(o), _s(_o.exceptions()) {
			_o.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit);
		}
		
		~ForceExceptions() {
			_o.exceptions(_s);
		}
	};

public:	
	PairwiseOrthogonalTransform(size_t x, size_t y, size_t z, bool bigEndian, bool isBSQ, bool isLossless)
		: _x(x), _y(y), _z(z), _bigEndian(bigEndian), _isBSQ(isBSQ), _isLossless(isLossless) { }
	
	void apply (istream & inputStream, ostream & outputStream, ostream & sideInformationStream) {
		// Setup input/output operations
		ForceExceptions fe1(inputStream);
		ForceExceptions fe2(outputStream);
		ForceExceptions fe3(sideInformationStream);
		
		FileOrder<FileFormatType> * in;
		FileOrder<FormatType> * out;
		
		if (_isBSQ) {
			in = new BSQ<FileFormatType>(_z, _y, _x);
			out = new BSQ<FormatType>(_z, _y, _x);
		} else {
			in = new BIL<FileFormatType>(_z, _y, _x);
			out = new BIL<FormatType>(_z, _y, _x);
		}
				
		BandLineReader<FormatType, FileFormatType> input(inputStream, _x, *in, _bigEndian);
		BandLineWriter<FormatType, FormatType> output(outputStream, _x, *out, _bigEndian);
		
		// Create the structure
		createStructure(&output);
		
		int firstLevelStart = pairs.size() - (_z / 2);
		
		// Read bands with BandReader and flush them into the structure
		for (size_t row = 0; row < _y; row++) {
			
			for (size_t band = 0; band < _z; band++) {
				FormatType * bl = input.receive(band);

				applyZeroMean (bl);

				if (band != _z - 1 || _z % 2 == 0) {
					pairs[firstLevelStart + band / 2]->send(bl, band % 2);
				} else {					
					/* this corresponds with the last pair of the previous level,
					 * since in the second level unpaired components are aligned left. */
					if (_z > 1) {
						// If there is more than one level
						pairs[firstLevelStart - 1]->send(bl, 1);
					} else {
						// Otherwise flush it to the output directly
						output.send(bl, 0);
					}
				}
			}
			
			input.startNextLine();
			output.startNextLine();
			sideInformation.writeToStream(sideInformationStream);
		}
		
		// Clean up
		destroyStructure();
		
		delete in;
		delete out;
	}
	
	void remove (istream & inputStream, ostream & outputStream, istream & sideInformationStream) {
		// Setup input/output operations
		ForceExceptions fe1(inputStream);
		ForceExceptions fe2(outputStream);
		ForceExceptions fe3(sideInformationStream);
		
		FileOrder<FormatType> * in;
		FileOrder<FileFormatType> * out;
		
		if (_isBSQ) {
			in = new BSQ<FormatType>(_z, _y, _x);
			out = new BSQ<FileFormatType>(_z, _y, _x);
		} else {
			in = new BIL<FormatType>(_z, _y, _x);
			out = new BIL<FileFormatType>(_z, _y, _x);
		}
		
		BandLineReader<FormatType, FormatType> input(inputStream, _x, *in, _bigEndian);
		BandLineWriter<FormatType, FileFormatType> output(outputStream, _x, *out, _bigEndian);
		
		// Create the structure
		createStructure(&input);
		
		int firstLevelStart = pairs.size() - (_z / 2);
		
		// Read bands from the structure and flush them into a BandWriter
		for (size_t row = 0; row < _y; row++) {
			sideInformation.readFromStream(sideInformationStream, _z);
			
			for (size_t band2 = 0; band2 < _z; band2++) {
				FormatType * bl;
				size_t band = _z - band2 - 1; 
				
				if (band != _z - 1 || _z % 2 == 0) {
					bl = pairs[firstLevelStart + band / 2]->receive(band % 2);
				} else {
					if (_z > 1) {
						bl = pairs[firstLevelStart - 1]->receive(1);
					} else {
						bl = input.receive(0);
					}
				}
	
				removeZeroMean(bl);			
				output.send(bl, band);
			}
			
			input.startNextLine();
			output.startNextLine();
		}
		
		// Clean up
		destroyStructure();
		
		delete in;
		delete out;
	}
}; 

}

#endif // POT_H 
