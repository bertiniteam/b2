
/** Define hardware precision to be 32 or 64 bits (or more?). Dependent on system and easy enough to change. */
#ifndef HARDWARE_PRECISION
#define HARDWARE_PRECISION 64
#endif

#ifndef __amp_hpp__
#define __amp_hpp__
#include <gmp.h>
#include <mpfr.h>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "Symbol.h"




class AMP : public Symbol
{
public:

  /** Constructors */

  /** Default Constructor.  Allocate all pointers, precision, and initalize the dereferenced mpfr_t pointer. */
  AMP();
  /** Double constructor.
      @param double v - The value of the dereferenced pointer. 
  */
  AMP(double v);
  
  /** MP String Constructor 
   @param std::string - the string of the constructor.
   @param int base - the base of the string
   @param prec - the precision in bits
   */
  AMP(std::string str, int base, int prec);
  /** MP mpfr_t constructor
      @param mpfr_t m - The new value of the derefenced mp.
  */
  AMP(mpfr_t m);


  /** The default destructor. */
  ~AMP();

  /** Copy constructors. */

  /** Reference copy constructor. */
  AMP(AMP& copy_from_me);
  /** Constant reference copy constructor. */
  AMP(const AMP& copy_from_me);
  /** Volatile reference copy constructor. */
  AMP(volatile AMP& copy_from_me);
  /** Constant volatile reference copy constructor. */
  AMP(const volatile AMP& copy_from_me);

  // getters

  /** Get the precision.
      @return precision - The precision in bits. 
   */
  int GetPrecision(); 

  /** Get the pointer to the double.
      @return d - The pointer to the double. 
  */
  inline double* GetD(){ return d; }

  /** Get the pointer to the mpfr object. 
      @return mp - The pointer to the mpfr object.
  */
  inline mpfr_t* GetMP(){ return mp; }

  /** A function for testing purposes that returns a string object with the AMP object's information.
      @return std::string - The string that contains the information.
  */
  std::string GetStats(std::string s);


  /** A function that returns true if using multiprecision and false when using hardware precision. */
  inline  bool GetPrecisionB(){
  //    std::cout << "int(mpfr_get_prec(*mp)) = ";
  //    std::cout << int(mpfr_get_prec(*mp)) << "\n";
      return (int(mpfr_get_prec(*mp) > HARDWARE_PRECISION ? true : false)); }

  // modifiers

  /** Set the precision. 
      @param - The precision in bits.
  */  
  void SetPrecision(int prec);


  /**
   Set the precision.        
   @param mpfr_prec mp_prec - The precision in bits as an mpfr_prec object.    
  */
  void SetPrecision(mp_prec_t mp_prec);

  /** Set the dereferenced mpfr_t pointer.  Uses the mpfr_set_str(mpfr_t m, char* c, int base, int prec) function from mpfr.
      @param std::string str - The floating represented as a string object. 
      @param int base - The base (i.e. number of symbols) of the string object
      @param prec - The precision in bits.
  */



  void Set(std::string str, int base, int prec);
  /** Set the dereferenced mpfr_t pointer. 
      @param mpfr_t m - The new value of the dereferenced pointer.
  */
  void Set(mpfr_t m);

  /** Set the dereferenced double pointer.
      @param double t - The value of the dereferenced pointer.
  */
  void Set(double t);
  
  /** The assignment = operator. 
   @param AMP a - The value to be assigned to *this. */
  AMP& operator = (AMP a);


  /** Addition operator. */

  friend AMP operator + (AMP a, AMP b);

  friend AMP operator + (AMP a, double b);
  
  friend AMP operator + (double a, AMP b);
  
  /** Subtraction operator. */

  friend AMP operator - (AMP a, AMP b);

  friend AMP operator - (AMP a, double b);

  friend AMP operator - (double a, AMP b);
  /** Multiplication operator. */

  friend AMP operator * (AMP a, AMP b);
  friend AMP operator * (double a, AMP b);
  friend AMP operator * (AMP a, double b);

  /** Division operator. */

  friend AMP operator / (AMP a, AMP b);
  friend AMP operator / (AMP a, double b);
  friend AMP operator / (double a, AMP b);

    
    
    /*  Inherited from Symbol   */
    virtual std::stringstream print(){std::stringstream ret; ret<< *(this->GetD()); return ret;};
    virtual Symbol* add(Symbol* operand);
    virtual Symbol* sub(Symbol* operand);
    virtual Symbol* mult(Symbol* operand);
    virtual Symbol* exp(int exp);
    virtual Symbol* neg();
    
    
    virtual AMP* clone ()
    {
        return new AMP(*this);
    }
    
  //protected:


private:
  /** Flag for using multiprecision. */ 
  //  bool precision;
  /** The pointer to the  mpfr_t object. */
  mpfr_t* mp;
  /** The pointer to the double. */
  double* d;
  
};


std::ostream& operator << (std::ostream &sout, mpfr_t m);

/** Insertion operator. */

std::ostream& operator << (std::ostream &sout, AMP a);

/** Extraction operator. */

//std::istream& operator >> (std::istream &sin, AMP a);




#endif
