
/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#ifndef _COMPLEMENTIF_H_
#define _COMPLEMENTIF_H_

#include "AMReX_RealVect.H"
#include "AMReX_BaseIF.H"

namespace amrex
{
    
  ///
  /**
     This implicit function is the negative of the implicit function with which
     it is constructed.  This results in an implicit function whose interior is
     the complement of the interior of the given implicit function (minus the
     the surface).
  */
  class ComplementIF: public BaseIF
  {
  public:
    ///
    /**
       Constructor specifying the implicit function and whether to complement
       (default = true)
    */
    ComplementIF(const BaseIF& a_impFunc);

    
    /// Destructor
    virtual ~ComplementIF();
    
    ///
    /**
       Return the value of the function at a_point.
    */
    virtual Real value(const RealVect& a_point) const;
    

    ///
    virtual BaseIF* newImplicitFunction() const;
    
    
  protected:

    
    BaseIF* m_impFunc; // implicit function to complement
    
  private:

    ComplementIF();
    void operator=(const ComplementIF& a_inputIF);
  };
}

#endif
