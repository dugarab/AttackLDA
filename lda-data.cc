// (C) Copyright 2004, David M. Blei (blei [at] cs [dot] cmu [dot] edu)

// This file is part of LDA-C.

// LDA-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// LDA-C is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#include "lda-data.h"

namespace attack_lda {
  void ReadDictionary( Corpus* corpus, char* dict_filename )
  {
    // statistics of the total number of words in dictionary
    FILE* fin = fopen( dict_filename , "r" );
    char tmp[300];
    corpus->num_terms = 0;
    while( fscanf( fin , "%s\n" , tmp) == 1 ) corpus->num_terms += 1;
    fclose(fin);

    // read the words in dictionary
    fin = fopen( dict_filename , "r" );
    corpus->dict = (char**) malloc(sizeof(char*) * corpus->num_terms);
    int v = 0;
    for(v = 0; v < corpus->num_terms; v++)
    {
      corpus->dict[v] = (char*) malloc(sizeof(char) * 300);
      fscanf( fin , "%s\n" , corpus->dict[v] );
    }
    fclose(fin);
  }


  Corpus* ReadCorpus(char* data_filename, char* dict_filename)
  {
    FILE *fileptr;
    int length, count, word, n, nd;
    Corpus* c;

    printf("reading data from %s\n", data_filename);
    c = new Corpus();
    c->docs = 0;
    ReadDictionary( c, dict_filename );
    c->num_docs = 0;
    fileptr = fopen(data_filename, "r");
    nd = 0;
    while ((fscanf(fileptr, "%10d", &length) != EOF))
    {
      c->docs = (Document*) realloc(c->docs, sizeof(Document)*(nd+1));
      c->docs[nd].length = length;
      c->docs[nd].total = 0.0;
      c->docs[nd].attack_total = 0.0;
      c->docs[nd].words = new std::vector<int>();
      c->docs[nd].words->resize( length );
      c->docs[nd].counts = new std::vector<double>();
      c->docs[nd].counts->resize( length );
      for (n = 0; n < length; n++)
      {
        fscanf(fileptr, "%10d:%10d", &word, &count);
        word = word - OFFSET;
        (*(c->docs[nd].words))[n] = word;
        (*(c->docs[nd].counts))[n] = count ;	
        c->docs[nd].total += count;
      }
      nd++;
    }
    fclose(fileptr);
    c->num_docs = nd;

    printf("number of docs    : %d\n", c->num_docs);
    printf("number of terms   : %d\n", c->num_terms);
    return(c);
  }

int MaxDocumentLength(Corpus* c)
{
  int n, max = 0;
  for (n = 0; n < c->num_docs; n++)
  if (c->docs[n].length > max) max = c->docs[n].length;
  return(max);
}
} // namespace attack_lda
