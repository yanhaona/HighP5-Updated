/* File: parser.h
 * --------------
 * This file provides constants and type definitions that will
 * are used and/or exported by the yacc-generated parser.
 */

#ifndef _H_parser
#define _H_parser

// here we need to include things needed for the yylval union
  
#include "scanner.h"            		// for MaxIdentLen
#include "../src/utils/list.h"
#include "../src/syntax/ast.h"
#include "../src/syntax/ast_type.h"
#include "../src/syntax/ast_def.h"
#include "../src/syntax/ast_task.h"
#include "../src/syntax/ast_partition.h"
#include "../src/syntax/ast_stmt.h"
#include "../src/syntax/ast_expr.h"
#include "../src/syntax/ast_library_fn.h"
#include "../src/semantics/task_space.h"		// for PartitionOrder
 
// Next, we want to get the exported defines for the token codes and
// typedef for YYSTYPE and exported global variable yylval.  These
// definitions are generated and written to the y.tab.h header file. But
// because that header does not have any protection against being
// re-included and those definitions are also present in the y.tab.c,
// we can get into trouble if we don't take precaution to not include if
// we are compiling y.tab.c, which we use the YYBISON symbol for. 
// Managing C headers can be such a mess! 

#ifndef YYBISON                 
#include "../y.tab.h"              
#endif

int yyparse();              // Defined in the generated y.tab.c file
void InitParser();          // Defined in parser.y

#endif
