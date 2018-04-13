# compressionStreams Notes

## Some key questions

* What C++ standard, 98, 11, 17?
* What about C printf/scanf implementation?
* Where to get STL lib source (primarily for review)?
* Probably any compression alg ==> binary stream
* What if used on one of the standard streams (cin, cout, cerr)?
  * doesn't make sense to do those 'cause they are ascii, right?

## Existing manipulators (do we need to also define funcs?)

### No parameters

* do all these work on << and >>?
* boolalpha ('1' vs 'True')
* showbase (C++ literal constant prefixes ('0x' for hex, '0' for octal)
* showpoint (always show decimal point)
* showpos (always show +/- sign)
* skipws (skip whitespace)
* unitbuf (flush upon completion of each insert operation)
* uppercase (all upper case)
* dec,hex,oct (format for integers)
* fixed/scientific (3.1415962 vs 3.1415E+000)
* internal/left/right (adjustment flags)
* ws (input manip. to extract whitespace)
* endl/ends/flush (output manipulators)

### Maniplators with parameters

* setiosflags(ios_bse::fmt_flags &flags) set all flags in one manipulator insert
* setbase(int base) (8,10,16)
* setfill(char f)
* setprecision(int n)
* setw(int width) field width

# Possible interface for compressor

* setzalag(char const \*algname)
  * sets the name of the compression algorithm
  * setzalag(0) ends the use of the alg for compression
* setzparams(char const \*)
  * sets compression algorithm params
  * same interface for all algorithms
  * each alg must parse out of the string the params

* insert operator won't necessarily present data to stream in order compresser needs it
* zfp wants 4^d blocks so either caller must present data in that order or stream needs
  to buffer it so that it can extract such order
  * buffer size/shape maybe settable as a param on the stream
