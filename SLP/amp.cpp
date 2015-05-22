#include <iomanip>
#include <string>
#include "amp.hpp"
#include <cmath>

//#define DEBUGGING 1

/** Constructors */

/** Default Constructor.  Allocate all pointers, precision, and initalize the dereferenced mpfr_t pointer. */


AMP::AMP() :  d(new double), mp(new mpfr_t[1]){
    //Set Symbol vars
    isVar = 0;
  
  *d = 0.0;
  mpfr_init2(*mp, HARDWARE_PRECISION);

  mpfr_set_str(*mp, "0.0", 10, MPFR_RNDN);


  
}
/** Double constructor.
    @param double v - The value of the dereferenced pointer. 
*/
AMP::AMP(double v): d(new double), mp(new mpfr_t[1]){
    //Set Symbol vars
    isVar = 0;
    
  *d = v;
  mpfr_init2(*mp, HARDWARE_PRECISION);
    std::stringstream ss;
    ss << std::setprecision(16)<< *d;
    
 
    mpfr_set_str(*mp, ss.str().c_str(), 10, MPFR_RNDN);
}
  
/** MP String Constructor 
    @param std::string - the string of the constructor.
    @param int base - the base of the string
    @param prec - the precision in bits
*/
AMP::AMP(std::string str, int base, int prec):d(new double), mp(new mpfr_t[1]){
    //Set Symbol vars
    isVar = 0;
    
#if DEBUGGING
    std::cout << "Entering AMP(std::string, int, int) Constructor\n";
#endif
    

  mpfr_init2(*mp, prec);
  mpfr_set_str(*mp, str.c_str(), base, MPFR_RNDN);
  *d = mpfr_get_d(*mp, MPFR_RNDN);

    
#if DEBUGGING
    std::cout << "Leaving AMP(std::string, int, int) Constructor\n";
#endif
}
/** MP mpfr_t constructor.
    @param mpfr_t m -- the value of the dereferenced mp
*/
AMP::AMP(mpfr_t m): d(new double), mp(new mpfr_t[1]){

  mpfr_init2(*mp, mpfr_get_prec(m));
  mpfr_set(*mp, m, MPFR_RNDN);
  *d = mpfr_get_d(*mp, MPFR_RNDN);

}


  /** The default destructor. */
AMP::~AMP(){
    
  delete d;  
  mpfr_clear(*mp);
  delete[] mp;

}

  /** Copy constructors. */

  /** Reference copy constructor. */
AMP::AMP(AMP& copy_from_me) : d(new double), mp(new mpfr_t[1]) {
    isVar = 0;
 // mp = new mpfr_t[ (mpfr_get_prec( *copy_from_me.mp ) % 310) + 1];
    
  *(this->d) = *(copy_from_me.d);  
  mpfr_init2(*(this->mp), copy_from_me.GetPrecision());
  mpfr_set(*(this->mp), *(copy_from_me.GetMP()), MPFR_RNDN);


}
/** Constant reference copy constructor. */
AMP::AMP(const AMP& copy_from_me):  d(new double), mp(new mpfr_t[1]) {
    isVar = 0;

  *(this->d) = *(copy_from_me.d);
  mpfr_init2(*(this->mp), mpfr_get_prec(*copy_from_me.mp) );
  mpfr_set(*(this->mp), *(copy_from_me.mp), MPFR_RNDN);
  

}
/** Volatile reference copy constructor. */
AMP::AMP(volatile AMP& copy_from_me): d(new double), mp(new mpfr_t[1]) {
    isVar = 0;
    
  *(this->d) = *(copy_from_me.d);
  mpfr_init2(*(this->mp), mpfr_get_prec(*copy_from_me.mp));
  mpfr_set(*(this->mp), *(copy_from_me.mp), MPFR_RNDN);

}
/** Constant volatile reference copy constructor. */
AMP::AMP(const volatile AMP& copy_from_me): d(new double), mp(new mpfr_t[1]) {
    isVar = 0;
    
  *(this->d) = *(copy_from_me.d);
  mpfr_init2(*(this->mp), mpfr_get_prec(*copy_from_me.mp));
  mpfr_set(*(this->mp), *(copy_from_me.mp), MPFR_RNDN);
  

}

// getters

/** Get the precision.
    @return precision - The precision in bits.
*/
int AMP::GetPrecision(){ 
   
  if (GetPrecisionB()){   
    return int(mpfr_get_prec(*mp));
  }
  else{
    return HARDWARE_PRECISION;
  }
  
}





// modifiers

/** Set the precision. 
    @param int prec - The precision in bits.
*/  
void AMP::SetPrecision(int prec){
  
    mpfr_t tmp;
    mpfr_init2(tmp, prec);
    
    // use a stringstream to get
    // the representation of the current
    // *mp object
    // and then set this to the tmp object.
    // this avoids having 'junk'
    // digits due to rounding errors and is more exact than the
    // mpfr_set_d function
    
    std::stringstream ss;
    ss << *mp;
    
    
    mpfr_set_str(tmp, ss.str().c_str(), 10, MPFR_RNDN);
    mpfr_clear(*mp);
    mpfr_init2(*mp, prec);
    ss.str("");
    ss << tmp;
    mpfr_set_str(*mp, ss.str().c_str(), 10, MPFR_RNDN);
    mpfr_clear(tmp);
    *d = mpfr_get_d(*mp, MPFR_RNDN);



}

/**
   Set the precision.
   @param mpfr_prec mp_prec - The precision in bits as an mpfr_prec object.
   @param bool b - This should be true if int(mp_prec) < HARDWARE_PRECISION
 */
void AMP::SetPrecision(mp_prec_t mp_prec){

 
 
    
    
    mpfr_t tmp;
    mpfr_init2(tmp, mp_prec);
    
    // use a stringstream to get
    // the representation of the current
    // *mp object
    // and then set this to the tmp object.
    // this avoids having 'junk'
    // digits due to rounding errors and is more exact than
    // the mpfr_set_d function
    
    std::stringstream ss;
    ss << *mp;
    
    
    mpfr_set_str(tmp, ss.str().c_str(), 10, MPFR_RNDN);
    mpfr_clear(*mp);
    mpfr_init2(*mp, mp_prec);
    ss.str("");
    ss << tmp;
    mpfr_set_str(*mp, ss.str().c_str(), 10, MPFR_RNDN);
    mpfr_clear(tmp);
    *d = mpfr_get_d(*mp, MPFR_RNDN);
    
    
    
}

/** Set the dereferenced mpfr_t pointer.  Uses the mpfr_set_str(mpfr_t m, char* c, int base, int prec) function from mpfr.
    @param std::string str - The floating represented as a string object. 
    @param int base - The base (i.e. number of symbols) of the string object
    @param prec - The precision in bits.
*/
void AMP::Set(std::string str, int base, int prec){

#if DEBUGGING
  std::cout << "Entering AMP::Set(std::string, int, int)\n";
#endif
  mpfr_clear(*mp);
  mpfr_init2(*mp, prec);  
  mpfr_set_str(*mp, str.c_str(), base, MPFR_RNDN);
  *d = mpfr_get_d(*mp, MPFR_RNDN);
#if DEBUGGING
  std::cout << "Leaving AMP::Set(std::string, int, int)\n";
#endif

}
/** Set the dereferenced mpfr_t pointer. 
    @param mpfr_t m - The new value of the dereferenced pointer.
*/
void AMP::Set(mpfr_t m){

#if DEBUGGING
  std::cout << "Entering AMP::Set(mpfr_t m)\n";
#endif


    
  // need to also set the precision properly

  mpfr_clear(*mp);
  mpfr_init2(*mp, mpfr_get_prec(m));
  mpfr_set(*mp, m, MPFR_RNDN);
  *d = mpfr_get_d(*mp, MPFR_RNDN);
   
    
#if DEBUGGING
  std::cout << "Leaving AMP::Set(mpfr_t m)\n";
#endif
}

/** Set the dereferenced double pointer.
    @param double t - The value of the dereferenced pointer.
*/
void AMP::Set(double t){
  
  *d = t;
  mpfr_clear(*mp);
  mpfr_init2(*mp, HARDWARE_PRECISION);
  std::stringstream ss;

  ss << std::setprecision(16)<< *d;
  mpfr_set_str(*mp, ss.str().c_str(), 10, MPFR_RNDN);
}

/** A function for testing purposes.
    returns a string with object information.
*/

std::string AMP::GetStats(std::string s){

  std::stringstream the_ss;
  
  the_ss << s << ".GetPrecisionB() = " << GetPrecisionB() << " at ";
  the_ss << (GetPrecisionB() ? GetPrecision() : HARDWARE_PRECISION) << " bits.\n";
  the_ss << "*" << s << " = " << *this << "\n";

  return the_ss.str();


}

/** The assignment = operator. 
    @param AMP a - The value to be assigned to *this. 
*/
AMP& AMP::operator = (AMP a){


  // avoid self-assignment

  if (this != &a){

    // these should never be null . . . 
    if (this->mp == NULL){
      (this->mp) = new mpfr_t[1];
    }
    if (this->d == NULL){
      (this->d) = new double;
    }
    

    
    if (a.GetPrecisionB()){
      mpfr_clear(*mp);
      mpfr_init2(*mp, a.GetPrecision());
      mpfr_set(*mp, *a.GetMP(), MPFR_RNDN);
      *d = *a.GetD();
    }
    else{
   
      *d = *a.GetD();
      mpfr_clear(*mp);
      mpfr_init2(*mp, HARDWARE_PRECISION);
      mpfr_set(*mp, *a.GetMP(), MPFR_RNDN);
    }
    
  }


  return *this;

}

/* An mpfr insertion operator. */

std::ostream& operator << (std::ostream &sout, mpfr_t m){

  if (m != NULL){
    

    // int mp_prec = int(mpfr_get_prec(m));
    // assume a base of 10 to output

    // d = log10(2^b);

    int buf_size = floor(int(mpfr_get_prec(m))*log10(2)) ;
    
    
    char* the_buf = new char[buf_size];
    std::stringstream somess;
    somess << "%." << buf_size << "RNf";
    mpfr_sprintf(the_buf, somess.str().c_str() , m);
    
    sout << the_buf;
    delete[] the_buf;
  }

  return sout;


}

/** Insertion operator. */

std::ostream& operator << (std::ostream &sout, AMP a){

  if (a.GetPrecisionB()){    
    sout << *a.GetMP();
  }
  else{
    sout << std::setprecision(16) << *a.GetD();
  }  
  return sout;
  
}


/** Extraction operator. */

// std::istream& operator >> (std::istream &sin, AMP a);

/** Addition operator. */

AMP operator + (AMP a, AMP b){

#if DEBUGGING
  std::cout << "Entering AMP::operator + (AMP a, AMP b)\n";
#endif
  AMP toRet;
  
  if (a.GetPrecisionB() || b.GetPrecisionB()){

    //    std::cout << "\n\nmultiprecision + \n \n";
    
    // one or the other has multiprecision
    // get the max
    int prec1 = a.GetPrecision();
    int prec2 = b.GetPrecision();

    if (prec1 > prec2){
      b.SetPrecision(prec1);
      toRet.SetPrecision(prec1);
   
    }
    else if (prec1 < prec2){
      a.SetPrecision(prec2);
      toRet.SetPrecision(prec2);
   
    }  
    else{      
      toRet.SetPrecision(prec1); // they are the same here, so either works     
    }

    mpfr_add(*toRet.GetMP(), *a.GetMP(), *b.GetMP(), MPFR_RNDN); 	
    *toRet.GetD() = mpfr_get_d(*toRet.GetMP(), MPFR_RNDN);
    
  }
  else{    
    
    *toRet.GetD() = *a.GetD() + *b.GetD();
    std::stringstream ss;
    ss << std::setprecision(16) << *toRet.GetD();
    mpfr_set_str(*toRet.GetMP(), ss.str().c_str(), 10, MPFR_RNDN);
  
    
  }
  



#if DEBUGGING
  std::cout << "Leaving AMP::operator + (AMP a, AMP b)\n";
#endif

  return toRet;
}



AMP operator + (AMP a, double b){

  AMP toRet;
  if (a.GetPrecisionB()){


    mpfr_t tmp;
    int the_prec = a.GetPrecision();
    mpfr_init2(tmp, mpfr_get_prec(*a.GetMP()));
  
    std::stringstream ss;
    ss << std::setprecision(16) << b;
      
    mpfr_set_str(tmp, ss.str().c_str(), 10, MPFR_RNDN);
      

    mpfr_add(tmp, tmp, *a.GetMP(), MPFR_RNDN);
    toRet.Set(tmp);
    mpfr_clear(tmp);


  }
  else{

    *toRet.GetD() = b + (*a.GetD());
    std::stringstream ss;
    ss << std::setprecision(16) << *toRet.GetD();
    mpfr_set_str(*toRet.GetMP(), ss.str().c_str(), 10, MPFR_RNDN);
    
  }
  return toRet;
  
}

AMP operator + (double a, AMP b){
    return b + a;
}

/** Subtraction operator. */

AMP operator - (AMP a, AMP b){

  AMP toRet;
  
  if (a.GetPrecisionB() || b.GetPrecisionB()){


    // one or the other has multiprecision
    // get the max
    int prec1 = a.GetPrecision();
    int prec2 = b.GetPrecision();
   

    if (prec1 > prec2){
      b.SetPrecision(prec1);
      toRet.SetPrecision(prec1);
   
    }
    else if (prec1 < prec2){
      a.SetPrecision(prec2);
      toRet.SetPrecision(prec2);
   
    }  
    else{      
      toRet.SetPrecision(prec1); // they are the same here, so either works     
   
    }
    mpfr_sub(*toRet.GetMP(), *a.GetMP(), *b.GetMP(), MPFR_RNDN); 	
    *toRet.GetD() = mpfr_get_d(*toRet.GetMP(), MPFR_RNDN);
  }
  else{    

    *toRet.GetD() = *a.GetD() - *b.GetD();
    std::stringstream ss;
    ss << std::setprecision(16) << *toRet.GetD();
    mpfr_set_str(*toRet.GetMP(), ss.str().c_str(), 10, MPFR_RNDN);
  }
  
  return toRet;

}


/** Multiplication operator. */

AMP operator * (AMP a, AMP b){


  AMP toRet;
  


  if (a.GetPrecisionB() || b.GetPrecisionB()){

    // std::cout << "\n\nmultiprecision * \n\n";
    // one or the other has multiprecision
    // get the max
    int prec1 = a.GetPrecision();
    int prec2 = b.GetPrecision();
   
    // int max = (prec1 > prec2 ? prec1 : prec2);
    if (prec1 > prec2){
      b.SetPrecision(prec1);
      toRet.SetPrecision(prec1);
   
    }
    else if (prec1 < prec2){
      a.SetPrecision(prec2);
      toRet.SetPrecision(prec2);
   
    }  
    else{      
      toRet.SetPrecision(prec1); // they are the same here, so either works     
   
    }
    mpfr_mul(*toRet.GetMP(), *a.GetMP(), *b.GetMP(), MPFR_RNDN); 	
    *toRet.GetD() = mpfr_get_d(*toRet.GetMP(), MPFR_RNDN);
  }
  else{    


    *toRet.GetD() = (*a.GetD()) * (*b.GetD());
    std::stringstream ss;
    ss << std::setprecision(16) << *toRet.GetD();
    mpfr_set_str(*toRet.GetMP(), ss.str().c_str(), 10, MPFR_RNDN);
  }
  
  return toRet;


}

// this avoid doings a double constructor for an AMP object if subracting a double - AMP if the AMP object precision is less than HARDWARE_PREC

AMP operator - (double a, AMP b){

  AMP toRet;
  
  if (b.GetPrecisionB()){
    mpfr_t tmp;
    int the_prec = b.GetPrecision();
    mpfr_init2(tmp, mpfr_get_prec(*b.GetMP()));

      
    std::stringstream ss;
    ss << std::setprecision(16) << a;
      
    mpfr_set_str(tmp, ss.str().c_str(), 10, MPFR_RNDN);
    
      
    mpfr_sub(tmp, tmp, *b.GetMP(), MPFR_RNDN);

      // we use Set(tmp) here so as to correctly have
      // the correct precision for *toRet.GetMP() -- default is HARDWARE_PRECISION
    toRet.Set(tmp);
    mpfr_clear(tmp);


  }
  else{

    *toRet.GetD() = a - (*b.GetD());
    std::stringstream ss;
    ss << std::setprecision(16) << *toRet.GetD();
    mpfr_set_str(*toRet.GetMP(), ss.str().c_str(),10, MPFR_RNDN);
  }
  return toRet;
    
    
}


AMP operator - (AMP a, double b){

  AMP toRet;
  
  if (a.GetPrecisionB()){
    mpfr_t tmp;
    int the_prec = a.GetPrecision();
    
    mpfr_init2(tmp, mpfr_get_prec(*a.GetMP()));

    std::stringstream ss;
    ss << std::setprecision(16) << b;
      
    mpfr_set_str(tmp, ss.str().c_str(), 10, MPFR_RNDN);
    mpfr_sub(tmp, *a.GetMP(), tmp, MPFR_RNDN);

    // we use Set(tmp) here so as to correctly have
    // the correct precision for *toRet.GetMP() -- default is HARDWARE_PRECISION
      
    toRet.Set(tmp);
    mpfr_clear(tmp);

  }
  else{

    *toRet.GetD() = (*a.GetD()) - b;
    std::stringstream ss;
    ss << std::setprecision(16) << *toRet.GetD();
    mpfr_set_str(*toRet.GetMP(), ss.str().c_str(), 10, MPFR_RNDN);
      
  }
  return toRet;
}





// this avoids doing a constructor for an AMP object if multiplying a double * AMP if the AMP object precision is less than HARDWARE_PREC

AMP operator * (double a, AMP b){
  // multiplying an AMP by a double, double is the first operand.  test b if multiprecision
  AMP toRet;

  if (b.GetPrecisionB()){
    
    mpfr_t tmp;
    int the_prec = b.GetPrecision();
    mpfr_init2(tmp, mpfr_get_prec(*b.GetMP()));
    std::stringstream ss;
    ss << std::setprecision(16) << a;
      
    mpfr_set_str(tmp, ss.str().c_str(),10, MPFR_RNDN);
    mpfr_mul(tmp, tmp, *b.GetMP(), MPFR_RNDN);

    // we use Set(tmp) here so as to correctly have
    // the correct precision for *toRet.GetMP() -- default is HARDWARE_PRECISION
    toRet.Set(tmp);
    mpfr_clear(tmp);
    
  }
  else{// precision is not used, so use hardware multiplication
    

    *toRet.GetD() = a*(*b.GetD());
    std::stringstream ss;
    ss << std::setprecision(16) << *toRet.GetD();
    mpfr_set_str(*toRet.GetMP(), ss.str().c_str(), 10, MPFR_RNDN);
  }
  return toRet;

}

// this avoids creating an AMP object from a double object when multiplying an AMP*double (similar to previous * operator, but swap these)

AMP operator * (AMP a, double b){
  return b*a; // * is commutative
}



/** Division operator. */

AMP operator / (AMP a, AMP b){
  AMP toRet;
#if DEBUGGING
    std::cout << "Entering AMP operator / (AMP a, AMP b) \n";
#endif
  if (a.GetPrecisionB() || b.GetPrecisionB()){

    //    std::cout << "\n\nmultiprecision / \n\n";
    // one or the other has multiprecision
    // get the max
    int prec1 = a.GetPrecision();
    int prec2 = b.GetPrecision();
   
    // int max = (prec1 > prec2 ? prec1 : prec2);
    if (prec1 > prec2){
      b.SetPrecision(prec1);
      toRet.SetPrecision(prec1);
   
    }
    else if (prec1 < prec2){
      a.SetPrecision(prec2);
      toRet.SetPrecision(prec2);
   
    }  
    else{      
      toRet.SetPrecision(prec1); // they are the same here, so either works     
   
    }
    mpfr_div(*toRet.GetMP(), *a.GetMP(), *b.GetMP(), MPFR_RNDN); 	
    *toRet.GetD() = mpfr_get_d(*toRet.GetMP(), MPFR_RNDN);
  }
  else{    

    *toRet.GetD() = (*a.GetD()) / (*b.GetD());
    std::stringstream ss;
    ss << std::setprecision(16) << *toRet.GetD();
    mpfr_set_str(*toRet.GetMP(), ss.str().c_str(), 10, MPFR_RNDN);
  }
  
    
#if DEBUGGING
    std::cout << "Leaving AMP operator / (AMP a, AMP b) \n";
#endif
    return toRet;


}



AMP operator / (double a, AMP b){

  AMP toRet;
  if (b.GetPrecisionB()){
    mpfr_t tmp;
    mpfr_inits2(mpfr_get_prec(*b.GetMP()), tmp,  mpfr_ptr(0));
    std::stringstream ss;
    ss << std::setprecision(16) << a;
      
    mpfr_set_str(tmp, ss.str().c_str(), 10, MPFR_RNDN);
      

    mpfr_div(tmp, tmp, *b.GetMP(), MPFR_RNDN);
    toRet.Set(tmp);
      mpfr_clear(tmp);
  }
  else{
    *toRet.GetD() = a / (*b.GetD());
    std::stringstream ss;
    ss << std::setprecision(16) << *toRet.GetD();
    mpfr_set_str(*toRet.GetMP(), ss.str().c_str(), 10, MPFR_RNDN);
  }


  return toRet;
}


AMP operator / (AMP a, double b){

  AMP toRet;

  if (a.GetPrecisionB()){
    mpfr_t tmp;
    mpfr_inits2(mpfr_get_prec(*a.GetMP()), tmp,  mpfr_ptr(0));
    
    std::stringstream ss;
    ss << std::setprecision(16) << b;
    mpfr_set_str(tmp, ss.str().c_str(), 10, MPFR_RNDN);
    mpfr_div(tmp, *a.GetMP(),tmp, MPFR_RNDN);

    toRet.Set(tmp);
    mpfr_clear(tmp);
  }
  else{
    *toRet.GetD() = (*a.GetD()) / b;
    std::stringstream ss;
    ss << std::setprecision(16) << *toRet.GetD();
    mpfr_set_str(*toRet.GetMP(), ss.str().c_str(), 10, MPFR_RNDN);
  }
  
  return toRet;

}











/*****************************************************************
                Inherited from Symbol
 *****************************************************************/
Symbol* AMP::add(Symbol* operand)
{
    AMP thisAMP = AMP(*this);
    
    if(typeid(*operand) == typeid(DoubleSymb))
    {
        DoubleSymb* tempop;
        tempop = (DoubleSymb*)operand;
        AMP* retop = new AMP(thisAMP + tempop->getValue());
        return retop;
    }
    else if(typeid(*operand) == typeid(IntSymb))
    {
        IntSymb* tempop = (IntSymb*)operand;
        AMP* retop = new AMP(thisAMP + tempop->getValue());
        return retop;
    }
    else if(typeid(*operand) == typeid(AMP))
    {
        AMP tempop = *((AMP*)operand);
        AMP* retop = new AMP(thisAMP + tempop);
        return retop;
    }
    else
    {
        Symbol* retop = 0;
        return retop;
    }
    
}





Symbol* AMP::sub(Symbol* operand)
{
    AMP thisAMP = AMP(*this);
    
    if(typeid(*operand) == typeid(DoubleSymb))
    {
        DoubleSymb* tempop;
        tempop = (DoubleSymb*)operand;
        AMP* retop = new AMP(thisAMP - tempop->getValue());
        return retop;
    }
    else if(typeid(*operand) == typeid(IntSymb))
    {
        IntSymb* tempop = (IntSymb*)operand;
        AMP* retop = new AMP(thisAMP - tempop->getValue());
        return retop;
    }
    else if(typeid(*operand) == typeid(AMP))
    {
        AMP tempop = *((AMP*)operand);
        AMP* retop = new AMP(thisAMP - tempop);
        return retop;
    }
    else
    {
        Symbol* retop = 0;
        return retop;
    }
    
}









Symbol* AMP::mult(Symbol* operand)
{
    AMP thisAMP = AMP(*this);
    
    if(typeid(*operand) == typeid(DoubleSymb))
    {
        DoubleSymb* tempop;
        tempop = (DoubleSymb*)operand;
        AMP* retop = new AMP(thisAMP * tempop->getValue());
        return retop;
    }
    else if(typeid(*operand) == typeid(IntSymb))
    {
        IntSymb* tempop = (IntSymb*)operand;
        AMP* retop = new AMP(thisAMP * tempop->getValue());
        return retop;
    }
    else if(typeid(*operand) == typeid(AMP))
    {
        AMP tempop = *((AMP*)operand);
        AMP* retop = new AMP(thisAMP * tempop);
        return retop;
    }
    else
    {
        Symbol* retop = 0;
        return retop;
    }
    
}





Symbol* AMP::exp(int exp)
{
    return this;
}



Symbol* AMP::neg()
{
    AMP thisAMP = AMP(*this);
    
    return new AMP((-1.0)*thisAMP);
}







