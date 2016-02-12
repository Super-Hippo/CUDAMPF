/**/
#ifndef NVRTCOPTIONS
#define NVRTCOPTIONS

#include <string>
#include <sstream>
using namespace std;

static const string __ARCH__("-arch=");
static const string __RDC__("-rdc=");
static const string __DBUG__("-G");
static const string __LN__("-lineinfo");
static const string __MAXREG__("-maxrregcount=");
static const string __INCLUDE__("-include=");
static const string __DEFINE__("-D");


/* Self-defined */
static const string __RIB__("-DRIB=");
static const string __SIZE__("-DSIZE=");
static const string __TOTAL__("-DTOTAL=");
static const string __Q__("-DQ=");

/* Operator */
static const string __EQUAL__("=");
static const string __F__("false");
static const string __T__("ture");

/* Compute Capacity */
static const string __CC35__("compute_35");
static const string __CC30__("compute_30");
static const string __CC20__("compute_20");


#endif /* NVRTCOPTIONS */

/* concatenation */
static inline string get_option(string A, string B)
{	
	string ap;
	ap.append(A);
	ap.append(B);
	return ap;
}

/* overload */
static inline string get_option(string A)
{
	string ap;
	ap.append(A);
	return ap;
}

/* other types to string */
static inline string int2str(int int_value)
{
	string string_temp;
	stringstream ss;
	ss << int_value;
	ss >> string_temp;
	return string_temp;
}

/* overload */
static inline string int2str(unsigned int int_value)
{
	string string_temp;
	stringstream ss;
	ss << int_value;
	ss >> string_temp;
	return string_temp;
}
