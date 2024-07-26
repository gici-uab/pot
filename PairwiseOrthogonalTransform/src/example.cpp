// --- this file has a column width of 120 characters ---
/**
 * @file
 * Example code of the usage of the Pairwise Orthogonal Transform.
 * See the file @a pot.h for details. 
 */
 
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<map>

#include "pot.h"

using namespace std;

string basename(const string & name) {
	size_t s = name.rfind('/');
	if (s == string::npos)
		return name;
	else 
		return name.substr(s + 1);
}

/**
 * A short and dependency-less portable argument parser.
 */
class NaiveParser {
private:
	const int _argc;
	const char * const * _argv;

	int x, y, z;
	int format;
	bool foward;
	bool help;
	bool bsq;
	string inputFile, outputFile, sideInformationFile;

	bool popOne(int & i, int & value) {
		i++;
		
		if (i >= _argc)
			return false;
			
		istringstream iss(_argv[i]);
		
		return !(iss >> std::dec >> value).fail();
	}
	
public:

	// Getters
	int getZ() const { return z; }
	int getY() const { return y; }
	int getX() const { return x; }
	int getFormat() const { return format; }
	int isFoward() const { return foward; }
	int isHelp() const { return help; }
	int isBsq() const { return bsq; }
	const char * getInputFile() const { return inputFile.c_str(); }
	const char * getOutputFile() const { return outputFile.c_str(); }
	const char * getSideInformationFile() const { return sideInformationFile.c_str(); }

	NaiveParser(const int argc, const char * const * argv) : _argc(argc), _argv(argv),
		x(0), y(0), z(0), format(0), foward(true), help(false), bsq(false) { }

	bool parse() {
		vector<string> remainingOptions;
		bool moreOptions = true;
		
		for (int i = 1; i < _argc; i++) {
			if (moreOptions) {
				string option(_argv[i]);
				if (option.compare("--") == 0) {
					moreOptions = false;
				} else if (option.compare("-s") * option.compare("--size") == 0) {
					if (! popOne(i, z) || ! popOne(i, y) || ! popOne(i, x))
						return false;
				} else if (option.compare("-f") * option.compare("--format") == 0) {
					if (! popOne(i, format))
						return false;
						
					if (format < 0)
						return false;
				} else if (option.compare("-b") * option.compare("--bsq") == 0) {
					bsq = true;
				} else if (option.compare("-r") * option.compare("--remove") == 0) {
					foward = false;
				} else if (option.compare("-h") * option.compare("--help") == 0) {
					help = true;
					return true;
				} else {
					remainingOptions.push_back(_argv[i]);
				}
			} else {
				remainingOptions.push_back(_argv[i]);
			}
		}
		
		if (x <= 0 || y <= 0 || z <= 0)
			return false;
			
		if (remainingOptions.size() != 3)
			return false;
		
		inputFile = remainingOptions[0];
		outputFile = remainingOptions[1];
		sideInformationFile = remainingOptions[2];
		
		return true;
	}
	
	void showHelp() const {
		cout <<
		"Copyright (C) 2010 Ian Blanes.\n"
		"This is free software.  You may redistribute copies of it under the terms of\n"
		"the GNU Affero General Public License Version 3. There is NO WARRANTY, to the\n"
		"extent permitted by law.\n"
		"\n"
		"Usage: " << basename(_argv[0]) << " [options] input-file output-file side-information-file\n"
		"\n"
		"  -h, --help		print this message\n"
		"  -s, --size z y x	sets image size (mandatory)\n"
		"  -f, --format format	sets raw data format\n"
		"			  0 = s16be -> f32be (default)\n"
		"			  1 = s16be -> s32be\n"
		"			  2 = s16be -> f32le\n"
		"			  3 = s32be -> f64be\n"
		"			  4 = s32be -> s32be\n"
		"			  5 = u16be -> f32be\n"
		"			  6 = u16be -> s32be\n"
		"  -b, --bsq		reads data in bsq format (default is bil)\n"
		"  -r, --remove		removes the transform (default is to apply it)\n"
		"\n"
		"Reference: I. Blanes and J. Serra-Sagrista, \"Pairwise Orthogonal Transform for\n"
		"Spectral Image Coding,\" IEEE Trans. Geosci. Remote Sens., 2010, in press.\n"
		"\n"
		"Example: \n"
		"  " << basename(_argv[0]) << " -s 224 2206 614 cuprite-bil.raw pot.raw si.raw\n"
		"  " << basename(_argv[0]) << " -r -s 224 2206 614 pot.raw recovered.raw si.raw\n"
		"\n";
	}
};

ios_base::openmode in = ifstream::in | ifstream::binary;
ios_base::openmode out = ofstream::out | ofstream::binary | ofstream::trunc;

template <class FormatType, class FileFormatType, class MeanType, class RotationType>
static bool apply (const NaiveParser & naiveParser, bool bigEndian, bool isLossless) {
	ifstream inputStream(naiveParser.getInputFile(), in);
	ofstream outputStream(naiveParser.getOutputFile(), out);
	
	if (! inputStream.good() || ! outputStream.good()) {
		cerr << "File opening failed." << endl;
		return false;
	}
	
	if (naiveParser.isFoward()) {
		ofstream sideInformationStream(naiveParser.getSideInformationFile(), out);
	
		if (! sideInformationStream.good()) {
			cerr << "File opening failed." << endl;
			return false;
		}
	
		pot::PairwiseOrthogonalTransform<FormatType, FileFormatType, MeanType, RotationType> pot(naiveParser.getX(),
			naiveParser.getY(), naiveParser.getZ(), bigEndian, naiveParser.isBsq(), isLossless);
		pot.apply (inputStream, outputStream, sideInformationStream);
	
		sideInformationStream.close();		
	} else {
		ifstream sideInformationStream(naiveParser.getSideInformationFile(), in);

		if (! sideInformationStream.good()) {
			cerr << "File opening failed" << endl;
			return false;
		}
	
		pot::PairwiseOrthogonalTransform<FormatType, FileFormatType, MeanType, RotationType> pot(naiveParser.getX(),
			naiveParser.getY(), naiveParser.getZ(), bigEndian, naiveParser.isBsq(), isLossless);
		pot.remove (inputStream, outputStream, sideInformationStream);
	
		sideInformationStream.close();
	}

	inputStream.close();
	outputStream.close();
	
	return true;
}

int main (int argc, char *argv[]) {
	cout << "POT " << __DATE__ << " " << __TIME__ << endl;

	NaiveParser naiveParser(argc, argv);
	
	if (! naiveParser.parse()) {
		cout << "Command line parsing failed." << endl;
		cout << "Try `" << basename(argv[0]) << " --help' for more information." << endl;
		return 1;
	}
	
	if (naiveParser.isHelp()) {
		naiveParser.showHelp();
		return 0;
	}

	switch (naiveParser.getFormat()) {
		case 0: // s16be -> f32be
			//template <typename FormatType, typename FileFormatType, typename MeanType, typename RotationType>
			apply<float, short, short, unsigned short>(naiveParser, true /* big endian */, false /* lossy */);
			break;
		case 1: // s16be -> s32be
			apply<int, short, short, unsigned short>(naiveParser, true /* big endian */, true /* lossless */);
			break;
		case 2: // s16be -> f32le
			//template <typename FormatType, typename FileFormatType, typename MeanType, typename RotationType>
			apply<float, short, short, unsigned short>(naiveParser, false /* little endian */, false /* lossy */);
			break;
		case 3: // s32be -> f64be
			apply<double, int, int, unsigned short>(naiveParser, true /* big endian */, false /* lossy */);
			break;
		case 4: // s32be -> s32be
			apply<int, int, int, unsigned short>(naiveParser, true /* big endian */, true /* lossless */);
			break;
		case 5: // u16be -> f32be
			apply<float, unsigned short, int, unsigned short>(naiveParser, true /* big endian */, false /* lossy */);
			break;
		case 6: // u16be -> s32be
			apply<int, unsigned short, int, unsigned short>(naiveParser, true /* big endian */, true /* lossless */);
			break;
		//case 7: // your option here
		default:
			cout << "Invalid format selected." << endl;
			return 1;
	}
	
	return 0;
}
