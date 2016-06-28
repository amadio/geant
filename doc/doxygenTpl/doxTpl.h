#ifndef doxTpl_h
#define doxTpl_h

/**
   @file doxTpl.h Header file for doxTpl

   Usually there is not much to say about the file itself as it is better
   to link the information to the module or to the class
*/

/**
   @defgroup doxygenExample Doxygen example for geant

   This is the example group for the Doxygen example for geant
*/

/**
   @class     doxTpl
   @ingroup   doxygenExample
   @brief     Short description of the class
   @date      23-Jun-2016
   @author    John Doe
   @author    Jane Doe
   @version   1.0
   @note      This is just an example
   @warning   The code contains a bug
   @bug       This is the bug
   @todo      Fix it

   This is a somewhat longish description of the class. If you want to double
   a number here how you should do

~~~ {.cpp}
   int b=2;
   doxTpl tpl;
   int twob = oneMethod(b);
~~~

Here is an example of the result

| Input | Expected | Obtained|
| :---- | :----:   | ----:   |
| 1     | 2        | 2       |
| 2     | 4        | 4       |

A more accurate algorithm could be derived from the following formula

\f[
  |I_2|=\left| \int_{0}^T \psi(t)
   \left\{
      u(a,t)-
      \int_{\gamma(t)}^a
      \frac{d\theta}{k(\theta,t)}
      \int_{a}^\theta c(\xi)u_t(\xi,t)\,d\xi
   \right\} dt
\right|
\f]

For any further information you can ask our expert.

@image html assistant.jpg "Elegant Ruby de la Grande Charriere"
@image latex assistant.eps "Elegant Ruby de la Grande Charriere" width=4cm


*/
class doxTpl {
public:
  /**
     @brief This is the default constructor

     Be careful not to allocate any space in the default constructor if you
     are using ROOT, because this will lead to a leak.
  */
  doxTpl();

  /**
     @brief Duplicate input number.
     @todo  Make it constant

     This method duplicates the input number.
  */
  int oneMethod(int a);

private:
  int fMember1;    ///<! ROOT transient member
  int *fMember2;   ///<[fMember1] ROOT array
  double fMember3; ///< Another class member
};

#endif
