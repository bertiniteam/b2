/***

My_Class first;           // initialization by default constructor
My_Class second(first);   // initialization by copy constructor
My_Class third = first;   // Also initialization by copy constructor
second = third;           // assignment by copy assignment operator

*/



/** Define hardware precision to be 32 or 64 bits (or more?). Dependent on system and easy enough to change. */
//#ifndef HARDWARE_PRECISION
//#define HARDWARE_PRECISION 64
//#endif

//#ifndef __complex_amp_hpp__
//#define __complex_amp_hpp__

// #define DEBUGGING 1

#include <gmp.h>
#include <mpfr.h>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "complex_amp.hpp"





//class ComplexAMP{
//public:

  /** Constructors */

  /** Default Constructor.  Allocate all pointers. */
ComplexAMP::ComplexAMP(): real_amp(new AMP), imag_amp(new AMP){

#ifdef DEBUGGING
  std::cout << "Entering and Leaving Default ComplexAMP Constructor.\n";
#endif
  
}


  /** Double constructor.
      @param AMP r - The value of the dereferenced AMP* real_amp.
      @param AMP i - The value of the dereferenced AMP* imag_amp.
  */
ComplexAMP::ComplexAMP(AMP r, AMP i): real_amp(new AMP), imag_amp(new AMP){

#ifdef DEBUGGING
  std::cout << "Entering ComplexAMP(AMP r, AMP i)  Constructor.\n";
#endif

  *real_amp = r;
  *imag_amp = i;


#if DEBUGGING
  std::cout << "Leaving ComplexAMP(AMP r, AMP i)  Constructor.\n";
#endif

}
  
/** MP mpfr_t constructor.
    @param mpfr_t mp_r - The new value of the dereferenced AMP* real_amp.
    @param mpfr_t mp_i - The new value of the dereferenced AMP* imag_amp.
*/
ComplexAMP::ComplexAMP(mpfr_t mp_r, mpfr_t mp_i) : real_amp(new AMP), imag_amp(new AMP){
  
#if DEBUGGING
  std::cout << "Entering ComplexAMP(mpfr_t mp_r, mpfr_t mp_i)  Constructor.\n";
#endif

  real_amp->Set(mp_r);
  imag_amp->Set(mp_i);


#if DEBUGGING
  std::cout << "Leaving ComplexAMP(mpfr_t mp_r, mpfr_t mp_i)  Constructor.\n";
#endif

}

/** double constructor. 
    @param double r - The new value of the dereferenced AMP* real_amp.
    @param double i - The new value of the dereferenced AMP* imag_amp.
*/

ComplexAMP::ComplexAMP(double r, double i) : real_amp(new AMP), imag_amp(new AMP){
  

#if DEBUGGING
  std::cout << "Entering ComplexAMP(double r, double i)  Constructor.\n";
#endif

  real_amp->Set(r);
  imag_amp->Set(i);


#if DEBUGGING
  std::cout << "Leaving ComplexAMP(double r, double i)  Constructor.\n";
#endif

}


/** The default destructor. */
ComplexAMP::~ComplexAMP(){

#if DEBUGGING
  std::cout << "Entering ~ComplexAMP Destructor.\n";
#endif

  if (real_amp != NULL){
    delete real_amp; real_amp = NULL;
  }
  if (imag_amp != NULL){
    delete imag_amp; imag_amp = NULL;
  }

#ifdef DEBUGGING
  std::cout << "Leaving ~ComplexAMP Destructor.\n";
#endif

}

  /** Copy constructors. */

 /** Reference copy constructor. */
ComplexAMP::ComplexAMP(ComplexAMP& copy_from_me){
  
#ifdef DEBUGGING
  std::cout << "Entering ComplexAMP(ComplexAMP& copy_from_me) Copy Constructor.\n";
#endif


  
  // is this a possible issue with memory leaks ... ?                                                                                                                                                    
  (this->real_amp) = new AMP;
  (this->imag_amp) = new AMP;
  *(this->real_amp) = *(copy_from_me.real_amp);
  *(this->imag_amp) = *(copy_from_me.imag_amp);
    

#if DEBUGGING
  std::cout << "Leaving ComplexAMP(ComplexAMP& copy_from_me) Copy Constructor.\n";
#endif
 
  

}

  /** Constant reference copy constructor. */
ComplexAMP::ComplexAMP(const ComplexAMP& copy_from_me){
  // is this a possible issue with memory leaks ... ?                                                                                                                                                  

#if DEBUGGING
  std::cout << "Entering ComplexAMP(const ComplexAMP& copy_from_me) Copy Constructor.\n";
#endif

  
  (this->real_amp) = new AMP;
  (this->imag_amp) = new AMP;
  *(this->real_amp) = *(copy_from_me.real_amp);
  *(this->imag_amp) = *(copy_from_me.imag_amp);

#if DEBUGGING
  std::cout << "Leaving ComplexAMP(const ComplexAMP& copy_from_me) Copy Constructor.\n";
#endif

}

/** Volatile reference copy constructor. */
ComplexAMP::ComplexAMP(volatile ComplexAMP& copy_from_me){
  // is this a possible issue with memory leaks ... ?                                                                                                                                                 

#if DEBUGGING
  std::cout << "Entering ComplexAMP(volatile ComplexAMP& copy_from_me) Copy Constructor.\n";
#endif

   
  (this->real_amp) = new AMP;
  (this->imag_amp) = new AMP;
  *(this->real_amp) = *(copy_from_me.real_amp);
  *(this->imag_amp) = *(copy_from_me.imag_amp);

#if DEBUGGING
  std::cout << "Leaving ComplexAMP(volatile ComplexAMP& copy_from_me) Copy Constructor.\n";
#endif

}

/** Constant volatile reference copy constructor. */
ComplexAMP::ComplexAMP(const volatile ComplexAMP& copy_from_me){
#if DEBUGGING
  std::cout << "Entering ComplexAMP(const volatile ComplexAMP& copy_from_me) Copy Constructor.\n";
#endif


  // is this a possible issue with memory leaks ... ?                                                                                                                                                    
  (this->real_amp) = new AMP;
  (this->imag_amp) = new AMP;
  *(this->real_amp) = *(copy_from_me.real_amp);
  *(this->imag_amp) = *(copy_from_me.imag_amp);

#if DEBUGGING
  std::cout << "Leaving ComplexAMP(const volatile ComplexAMP& copy_from_me) Copy Constructor.\n";
#endif

}

AMP ComplexAMP::NormSquared(){
  
 
  AMP toRet = *(this->GetReal())*(*(this->GetReal())) + *(this->GetImag())*(*(this->GetImag()));


  return toRet;
}

ComplexAMP ComplexAMP::Conjugate(){


  
  ComplexAMP toRet(*this); // copy this
  
  toRet.SetImag( 0.0 - *(toRet.GetImag())); // this should work ... note that 0.0 - *(toRet.GetImag()) is better than *(toRet.GetImag()) * -1
  
  return toRet;
}

std::string ComplexAMP::GetStats(std::string s){
  
#if DEBUGGING
  std::cout << "Entering ComplexAMP::GetStats(std::string s)\n";
#endif 
  std::stringstream the_ss;
  the_ss << s << ".GetPrecisionRealB() = " << GetPrecisionRealB() << " at "; 
  the_ss << (GetPrecisionRealB() ? GetRealPrecision() : HARDWARE_PRECISION) << " bits.\n";
  the_ss << s << ".GetPrecisionImagB() = " << GetPrecisionImagB() << " at ";
  the_ss << (GetPrecisionImagB() ? GetImagPrecision() : HARDWARE_PRECISION) << " bits.\n";
  the_ss << "*" << s << ".GetReal() = " << *GetReal() << "\n";
  the_ss << "*" << s << ".GetImag() = " << *GetImag();
  
  //  std::cout << "the_ss = " << the_ss.str();
#if DEBUGGING
  std::cout << "Leaving ComplexAMP::GetStats(std::string s)\n";
#endif 

  return the_ss.str();

}


// modifiers

/** The assignment = operator. 
    @param AMP a - The value to be assigned to *this. 
*/
ComplexAMP& ComplexAMP::operator = (ComplexAMP a){
   
#ifdef DEBUGGING
  std::cout << "Entering ComplexAMP = assignment operator.\n";
#endif

  if (this != &a){
 


    // technically, these if statements should never be true....
    if (this->real_amp == NULL){
      std::cout << "The real_amp pointer is null -- in ComplexAMP = assignment operator --------- super suspicious !\n";
      this->real_amp = new AMP;
    }
    if (this->imag_amp == NULL){
      std::cout << "The imag_amp pointer is null --- in ComplexAMP = assignment operator -------- super suspicious !\n";
      this->imag_amp = new AMP;
    }

    *(this->real_amp) = *(a.real_amp);
    *(this->imag_amp) = *(a.imag_amp);
  
  }

#if DEBUGGING
  std::cout << "Leaving ComplexAMP = assignment operator.\n";
#endif

  return *this;
}


/** Addition operator. */

ComplexAMP operator + (ComplexAMP a, ComplexAMP b){

  ComplexAMP toRet;

  
  // these can be set on two different threads
  toRet.SetReal(*a.GetReal() + *b.GetReal());
  toRet.SetImag(*a.GetImag() + *b.GetImag());

  return toRet;
}

ComplexAMP operator + (ComplexAMP a, AMP b){

#if DEBUGGING  
  std::cout << "Entering ComplexAMP operator + (ComplexAMP a, AMP b)\n";
#endif

  ComplexAMP toRet;

  toRet.SetReal(*a.GetReal() + b);
  toRet.SetImag(*a.GetImag());
  
  

  
#if DEBUGGING
  std::cout << "Leaving ComplexAMP operator + (ComplexAMP a, AMP b)\n";
#endif
  return toRet;
}
ComplexAMP operator + (AMP a, ComplexAMP b){

  return b + a;
}

ComplexAMP operator + (ComplexAMP a, double b){
  
  ComplexAMP toRet;
  if (a.GetPrecisionRealB() || a.GetPrecisionImagB()){
    int arealprecision = a.GetRealPrecision();
    int aimagprecision = a.GetImagPrecision();
    int maxprec = (arealprecision >= aimagprecision ? arealprecision : aimagprecision);
    
    a.SetRealPrecision(maxprec);
    a.SetImagPrecision(maxprec);
    
    toRet.SetReal(*a.GetReal() + b);
    toRet.SetImag(*a.GetImag());
  }
  else{
    toRet.SetReal(*(*a.GetReal()).GetD() + b);
    toRet.SetImag(*a.GetImag());
  }
  return toRet;
}

ComplexAMP operator + (double a, ComplexAMP b){

  return b + a;
}

/** Subtraction operator. */

ComplexAMP operator - (ComplexAMP a, ComplexAMP b){
  ComplexAMP toRet;


  // these can be set on two different threads
  toRet.SetReal(*a.GetReal() - *b.GetReal());
  toRet.SetImag(*a.GetImag() - *b.GetImag());

  return toRet;
}

ComplexAMP operator - (ComplexAMP a, AMP b){
  
  ComplexAMP toRet;

  toRet.SetReal(*a.GetReal() - b);
  toRet.SetImag(*a.GetImag());
  return toRet;
}
ComplexAMP operator - (AMP a, ComplexAMP b){
  ComplexAMP toRet;
  toRet.SetReal(a - *b.GetReal());
  toRet.SetImag(0.0 - (*b.GetImag()));
  return toRet;
}

ComplexAMP operator - (ComplexAMP a, double b){
  ComplexAMP toRet;
  
  toRet.SetReal(*a.GetReal() - b);
  toRet.SetImag(*a.GetImag());

  return toRet;
}


ComplexAMP operator - (double a, ComplexAMP b){

  ComplexAMP toRet;
  toRet.SetReal(a - *b.GetReal());
  toRet.SetImag(0.0 - (*b.GetImag()));
  return toRet;
}


/** Multiplication operator. */

ComplexAMP operator * (ComplexAMP a, ComplexAMP b){

  //  k1 = c * (a + b)
  //  k2 = a * (d − c)
  //  k3 = b * (c + d)

  //  Real part = k1 − k3
  //  Imaginary part = k1 + k2.
  ComplexAMP toRet;


  // these three can be placed on different threads
  AMP k1 = *(b.GetReal())*( *(a.GetReal()) + *(a.GetImag()) );
  AMP k2 = *(a.GetReal())*( *(b.GetImag()) - *(b.GetReal()) );
  AMP k3 = *(a.GetImag())*( *(b.GetReal()) + *(b.GetImag()) );
  
  // need to sync threads here.
  
  // these two can be placed on different threads

  toRet.SetReal(k1 - k3);
  toRet.SetImag(k1 + k2);
  return toRet;
}


// complex number * a real number 
// don't need to do gaussian trick
// just multiply both components of complex # by the real #

ComplexAMP operator * (ComplexAMP a, AMP b){

  ComplexAMP toRet;
  if (a.GetPrecisionRealB() || a.GetPrecisionImagB() || b.GetPrecisionB()){
    // one of these has bigger than 64-bits of precision
    // get the max precision
    int arealprec = a.GetRealPrecision();
    int aimagprec = a.GetImagPrecision();
    int bprec = b.GetPrecision();
    
    int maxprec = arealprec;
    if (maxprec < aimagprec){
      maxprec = aimagprec;
    }
    if (maxprec < bprec){
      maxprec = bprec;
    }
    // set them all to the maxprecision
    a.SetRealPrecision(maxprec);
    a.SetImagPrecision(maxprec);
    b.SetPrecision(maxprec);
    toRet.SetReal( (*a.GetReal()) * b);
    toRet.SetImag( (*a.GetImag()) * b); 
  }   
  else{
    toRet.SetReal( *(*a.GetReal()).GetD() * (*b.GetD()));
    toRet.SetImag( *(*a.GetImag()).GetD() * (*b.GetD()));
  }

  
  return toRet;

}


ComplexAMP operator * (AMP a, ComplexAMP b){
  return b*a;
}


ComplexAMP operator *(double a, ComplexAMP b){

#if DEBUGGING
  std::cout << "Entering operator * (double a, ComplexAMP b)\n";
#endif


 ComplexAMP toRet;
 if (b.GetPrecisionRealB() || b.GetPrecisionImagB()){
   mpfr_t tmp, tmp1, tmp2;
   
   int brealprec = b.GetRealPrecision();
   int bimagprec = b.GetImagPrecision();
   int maxprec = (brealprec >= bimagprec ? brealprec : bimagprec);
   
   mpfr_inits2(maxprec, tmp, tmp1, tmp2, mpfr_ptr(0));
   std::stringstream ss;
   ss  << std::setprecision(16) << a;
     
   mpfr_set_str(tmp, ss.str().c_str(), 10, MPFR_RNDN);
   
   mpfr_mul(tmp1, tmp, *(*b.GetReal()).GetMP(), MPFR_RNDN);
   mpfr_mul(tmp2, tmp, *(*b.GetImag()).GetMP(), MPFR_RNDN);
   
   toRet.SetReal(tmp1);
   toRet.SetImag(tmp2);
   
   mpfr_clears(tmp, tmp1, tmp2, mpfr_ptr(0));
   
 }
 else{
   toRet.SetReal(a*( *(*b.GetReal()).GetD()));
   toRet.SetImag(a*( *(*b.GetImag()).GetD()));
 }

#if DEBUGGING
  std::cout << "Leavinging operator * (double a, ComplexAMP b)\n";
#endif
 
 return toRet;
 
}

ComplexAMP operator * (ComplexAMP a, double b){
  
#if DEBUGGING
  std::cout << "Entering and Leaving operator * (ComplexAMP a, double b)\n";
#endif

  return b*a; //multiplication is commutative!
}


/** Division operator. */

ComplexAMP operator / (ComplexAMP a, ComplexAMP b){


  // a = a + bi, b = c + di
  // then a/b = (a+bi)/(c+di) = (a+bi)(c-di) / (c^2 + d^2)

  // creating these two objects can be done on different threads

  ComplexAMP toRet = a*b.Conjugate();
  AMP NormSquared = b.NormSquared();

  

  // need to sync threads here



  // these two operations can be placed on different threads

  toRet.SetReal(*(toRet.GetReal()) / NormSquared);
  toRet.SetImag(*(toRet.GetImag()) / NormSquared);


  return toRet;
}

ComplexAMP operator / (ComplexAMP a, AMP b){

  ComplexAMP toRet;
  if (a.GetPrecisionRealB() || a.GetPrecisionImagB() || b.GetPrecisionB()){
    // we have multiprecision!
    int arealprecision = a.GetRealPrecision();
    int aimagprecision = a.GetImagPrecision();
    int bprecision = b.GetPrecision();
    int maxprec = arealprecision;
    if (maxprec < aimagprecision){
      maxprec = aimagprecision;
    }
    if (maxprec < bprecision){
      maxprec = bprecision;
    }
    a.SetRealPrecision(maxprec);
    a.SetImagPrecision(maxprec);
    b.SetPrecision(maxprec);
    
    toRet.SetReal(*a.GetReal() / b);
    toRet.SetImag(*a.GetImag() / b);

  }
  else{
    toRet.SetReal(*a.GetReal() / (*b.GetD()));
    toRet.SetImag(*a.GetImag() / (*b.GetD())); 
  }

  return toRet;

}



ComplexAMP operator / (AMP a, ComplexAMP b){

  // we don't save anything, since we are dividing a real AMP object by a ComplexAMP object
  ComplexAMP toRet = a*b.Conjugate();

  AMP NormSquared = b.NormSquared();
  
  toRet.SetReal( (*toRet.GetReal()) / NormSquared );
  toRet.SetImag( (*toRet.GetImag()) / NormSquared );

  return toRet;

}

ComplexAMP operator / (double a, ComplexAMP b){

  ComplexAMP toRet = a*b.Conjugate();
  AMP NormSquared = b.NormSquared();
  toRet.SetReal( (*toRet.GetReal()) / NormSquared);
  toRet.SetImag( (*toRet.GetImag()) / NormSquared);
  
  return toRet;
}

ComplexAMP operator / (ComplexAMP a, double b){
  ComplexAMP toRet;

  if (a.GetPrecisionRealB() || a.GetPrecisionImagB()){
    int arealprecision = a.GetRealPrecision();
    int aimagprecision = a.GetImagPrecision();
    int maxprec = (arealprecision >= aimagprecision ? arealprecision : aimagprecision);
    a.SetRealPrecision(maxprec);
    a.SetImagPrecision(maxprec);    
    
    toRet.SetReal(*a.GetReal() / b);
    toRet.SetImag(*a.GetImag() / b);
  }
  else{

    // these are faster as they are hardware double division (return doubles) and then using the SetReal(double) and SetImag(double) functions
    toRet.SetReal( *(*a.GetReal()).GetD() / b); 
    toRet.SetImag( *(*a.GetImag()).GetD() / b);
  }
  
  return toRet;
}


/** Insertion operator. */

std::ostream& operator << (std::ostream &sout, ComplexAMP a){


  sout << *a.GetReal() << " ";
  // see if we need to add a + or not

  if (   (* ( *a.GetImag() ).GetD()) >= 0){
    sout << " + ";
}


  sout << *a.GetImag() << "*I";
  return sout;
}

/** Extraction operator. */

//std::istream& operator >> (std::istream &sin, AMP a);





