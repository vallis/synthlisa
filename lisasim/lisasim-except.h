/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */

#ifndef _LISASIM_EXCEPT_H_
#define _LISASIM_EXCEPT_H_

// later can be made to inherit from a more generic class

class SynthLISAException {};

class ExceptionOutOfBounds : SynthLISAException {};
class ExceptionUndefined : SynthLISAException {};
class ExceptionWrongArguments : SynthLISAException {};
class ExceptionFileError : SynthLISAException {};
class ExceptionKeyboardInterrupt : SynthLISAException {};

#endif /* _LISASIM_EXCEPT_H_ */
