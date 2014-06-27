#ifndef ATTACKLDA_LDADATA_H_
#define ATTACKLDA_LDADATA_H_

#include <stdio.h>
#include <stdlib.h>

#include "lda.h"

#define OFFSET 0;                  // offset for reading data

namespace attack_lda {
Corpus* ReadCorpus(char* data_filename, char* dict_filename);
void ReadDictionary( Corpus* corpus, char* dict_filename );

int MaxDocumentLength(Corpus* c);
} // namespace attack_lda
#endif
