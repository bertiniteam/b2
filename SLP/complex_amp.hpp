
/** Define hardware precision to be 32 or 64 bits (or more?). Dependent on system and easy enough to change. */
//#ifndef HARDWARE_PRECISION
//#define HARDWARE_PRECISION 64
//#endif

#ifndef __complex_amp_hpp__
#define __complex_amp_hpp__


#include <gmp.h>
#include <mpfr.h>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "amp.hpp"



class ComplexAMP{
public:

  /** Constructors */
  
  /** Default Constructor.  Allocate all pointers, precision, and initalize the dereferenced mpfr_t pointer. */
  ComplexAMP();
  /** Double constructor.
      @param AMP r - The value of the dereferenced AMP* real_amp.
      @param AMP i - The value of the dereferenced AMP* imag_amp.
  */
  ComplexAMP(AMP r, AMP i);
  
  /** MP mpfr_t constructor.
      @param mpfr_t mp_r - The new value of the dereferenced AMP* real_amp.
      @param mpfr_t mp_i - The new value of the dereferenced AMP* imag_amp.
  */
  ComplexAMP(mpfr_t mp_r, mpfr_t mp_i);

  /** double constructor. 
      @param double r - The new value of the dereferenced AMP* real_amp.
      @param double i - The new value of the dereferenced AMP* imag_amp.
   */

  ComplexAMP(double r, double i);


  /** The default destructor. */
  ~ComplexAMP();

  /** Copy constructors. */

  /** Reference copy constructor. */
  ComplexAMP(ComplexAMP& copy_from_me);

  /** Constant reference copy constructor. */
  ComplexAMP(const ComplexAMP& copy_from_me);

  /** Volatile reference copy constructor. */
  ComplexAMP(volatile ComplexAMP& copy_from_me);

  /** Constant volatile reference copy constructor. */
  ComplexAMP(const volatile ComplexAMP& copy_from_me);

  ComplexAMP Conjugate();
  AMP NormSquared();

  // getters


  /** A function that returns a string object with all the information of the ComplexAMP object.
      Returns a string that stores the value of both the real_amp and imag_amp objects with their precision in bits.
      @param std::string s - the 'name' of the object used in the program.
      @return string object that contains the information regarding precision and value.
  */

  std::string GetStats(std::string s);

  /** Get the precision of the dereferenced AMP* r.
      @return precision - The precision in bits. 
   */
  inline int GetRealPrecision(){ return real_amp->GetPrecision(); } 

  /** Get the precision of the dereferenced AMP* i.
      @return precision - The precision in bits.
  */
  inline int GetImagPrecision(){ return imag_amp->GetPrecision(); }

  /** Get the pointer to the real AMP.
      @return d - The pointer to the real AMP. 
  */
  inline AMP* GetReal(){ return real_amp; }

  /** Get the pointer to the imaginary AMP. 
      @return mp - The pointer to the imaginary AMP.
  */
  inline AMP* GetImag(){ return imag_amp; }



  // modifiers

  /** Set the precision of the dereferenced AMP* real. 
      @param int prec - The precision in bits.
  */  
  void SetRealPrecision(int prec){real_amp->SetPrecision(prec);}
  
  /** Set the precision of the dereferenced AMP* imag.
      @param int prec - The precision in bits.
   */
  void SetImagPrecision(int prec){imag_amp->SetPrecision(prec);}

  /** Set the precision of both the dereferenced AMP* imag and AMP* real.
      @param int prec - The precision in bits.
  */
  void SetPrecision(int prec){ SetRealPrecision(prec); SetImagPrecision(prec); }

  /** Set the dereferenced AMP* real.  Uses the mpfr_set_str(mpfr_t m, char* c, int base, int prec) function from mpfr.
      @param std::string str - The floating represented as a string object. 
      @param int base - The base (i.e. number of symbols) of the string object
      @param prec - The precision in bits.
  */
  void SetReal(std::string str, int base, int prec){
    real_amp->Set(str, base, prec);  
  }


  /** Set the dereferenced AMP* real. 
      @param mpfr_t m - The new value of the dereferenced pointer.
  */
  void SetReal(mpfr_t m){
    real_amp->Set(m);
  }
  
  /** Set the dereferenced AMP* real to the double value.
      @param double t - The double value of the dereferenced AMP* real.
  */
  void SetReal(double t){real_amp->Set(t);}

  /** Set the dereferenced AMP* real to the AMP value.
      @param AMP a - The new value of the dereferenced AMP* real.
   */
  void SetReal(AMP a){
#if DEBUGGING
    std::cout << "Entering ComplexAMP::SetReal(AMP a)\n";
#endif 
   *real_amp = a;
 
#if DEBUGGING
    std::cout << "Leaving ComplexAMP::SetReal(AMP a)\n";
#endif  
  }

  /** Set the dereferenced AMP* imag.  Uses the mpfr_set_str(mpfr_t m, char* c, int base, int prec) function from mpfr.
      @param std::string str - The floating represented as a string object. 
      @param int base - The base (i.e. number of symbols) of the string object
      @param prec - The precision in bits.
  */
  void SetImag(std::string str, int base, int prec){imag_amp->Set(str, base, prec);}


  /** Set the dereferenced AMP* imag. 
      @param mpfr_t m - The new value of the dereferenced pointer.
  */
  void SetImag(mpfr_t m){imag_amp->Set(m);}


  /** Set the dereferenced AMP* imag to the double value.
      @param double t - The double value of the dereferenced AMP* imag.
  */
  void SetImag(double t){imag_amp->Set(t);}
  

  /** Set the dereferenced AMP* real to the AMP value.
      @param AMP a - The new value of the dereferenced AMP* real.
   */
  void SetImag(AMP a){

#if DEBUGGING
    std::cout << "Entering ComplexAMP::SetImag(AMP a)\n";
#endif
    *imag_amp = a;

#if DEBUGGING
    std::cout << "Leaving ComplexAMP::SetImag(AMP a)\n";
#endif 
 }

  /** The assignment = operator. 
   @param AMP a - The value to be assigned to *this. */
  ComplexAMP& operator = (ComplexAMP a);


  /** Addition operator. */

  friend ComplexAMP operator + (ComplexAMP a, ComplexAMP b);

  friend ComplexAMP operator + (ComplexAMP a, AMP b);
  friend ComplexAMP operator + (AMP a, ComplexAMP b);

  /** Subtraction operator. */

  friend ComplexAMP operator - (ComplexAMP a, ComplexAMP b);
  friend ComplexAMP operator - (ComplexAMP a, AMP b);
  friend ComplexAMP operator - (AMP a, ComplexAMP b);

  /** Multiplication operator. */

  friend ComplexAMP operator * (ComplexAMP a, ComplexAMP b);

  friend ComplexAMP operator * (AMP a, ComplexAMP b);
  friend ComplexAMP operator * (ComplexAMP a, AMP b);
  friend ComplexAMP operator * (double a, ComplexAMP b);
  friend ComplexAMP operator * (ComplexAMP a, double b);
  

  /** Division operator. */

  friend ComplexAMP operator / (ComplexAMP a, ComplexAMP b);

  friend ComplexAMP operator / (ComplexAMP a, AMP b);
  friend ComplexAMP operator / (AMP a, ComplexAMP b);
  friend ComplexAMP operator / (ComplexAMP a, double b);
  friend ComplexAMP operator / (double a, ComplexAMP b);
  

  //protected:
  inline  bool GetPrecisionRealB(){ return real_amp->GetPrecisionB(); }
  inline  bool GetPrecisionImagB(){ return imag_amp->GetPrecisionB(); }

private:

  AMP* real_amp;
  AMP* imag_amp;

  
};


std::ostream& operator << (std::ostream &sout, mpfr_t m);

/** Insertion operator. */

std::ostream& operator << (std::ostream &sout, ComplexAMP a);

/** Extraction operator. */

//std::istream& operator >> (std::istream &sin, AMP a);




#endif
