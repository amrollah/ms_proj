//
// Efficient CSV Reader for MATLAB
//
// Copyright (C) 2012, Stanislaw Adaszewski
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// minor modifications by Jan Poland, ABB.CH-RD.C1, 2015-02-18

#include <mex.h>
#include <vector>
#include <string>

#define NROWS_LIM 5000000

using namespace std;

static double zero = 0.0;
static double NaN = zero / 0;

static char sep[81];
static int nsep;
static char quote = '"';
static char escape = '\\';

void insert_elem(const char *nptr, int s_size, int m, int n, int M, double *p, mxArray *pcell) {
            char *endptr;
            double d = strtod(nptr, &endptr);
            if (s_size == 0 || endptr != nptr + s_size) {
                p[n * M + m] = NaN;
                if (pcell) {
                    mxArray *tmp = mxCreateString(nptr);
                    mxSetCell(pcell, (int)(n * M + m), tmp);
                }
            } else {
                p[n * M + m] = d;
            }
}    

void parse_file(FILE *f, vector<vector<string> > &data, size_t &max_columns, int &nrows, double *p, mxArray *pcell) {
    char buf[4096];
    int state = 0, isep, storing = 1, nrow = 0, ncol = 0;
    string cell;
    vector<string> row;
    if (p!=NULL) storing = 2;
    while (true) {
        size_t n = fread(buf, 1, 4096, f);
        if (n == 0) {
            break;
        }
        for (size_t i = 0; i < n; i++) {
            switch(state) {
                case 0: // virgin state
                    for (isep=0;isep<nsep;isep++) if (buf[i] == sep[isep]) break;
                    if (buf[i] == quote) {
                        state = 1; // quoted string
                    } else if (isep<nsep || buf[i] == '\n') { //in the separators or CR
                        if (storing==1) row.push_back(cell);
                        else if (storing==2) insert_elem(cell.c_str(),cell.size(),nrow,ncol,nrows,p,pcell);
                        cell = "";
                        ncol ++;
                    } else if (buf[i] == '\r') {
                        // do nothing
                    } else {
                        cell += buf[i];
                    }
                    if (buf[i] == '\n') {
                        if (ncol > max_columns) max_columns = ncol;
                        if (storing==1) {
                          if (nrow<NROWS_LIM) data.push_back(row);
                          else {storing = 0; data.clear();}
                          row.clear();
                        } else if (storing==2) {
                          for (; ncol < max_columns; ncol++) p[ncol * nrows + nrow] = NaN;
                        }
                        ncol = 0; nrow++;
                    }
                    break;
                case 1: // quoted string
                    if (buf[i] == escape) {
                        state = 2; // escaped state
                    } else if (buf[i] == quote) {
                        state = 3; // potential double quote
                    } else {
                        cell += buf[i];
                    }
                    break;
                case 2: // escaped
                    if (buf[i] == 'n') {
                        cell += '\n';
                    } else if (buf[i] == 't') {
                        cell += '\t';
                    } else if (buf[i] == 'r') {
                        cell += '\r';
                    } else if (buf[i] == quote) {
                        cell += quote;
                    } else {
                        mexErrMsgTxt("Encountered unknown escape sequence.");
                    }
                    state = 1;
                    break;
                case 3: // potential double quote
                    if (buf[i] == quote) {
                        cell += quote;
                        state = 1; // quoted string continues
                    } else {
                        state = 0; // parse it in the regular way
                        i--;
                    }
                    break;
            }
        }
    }
    nrows = nrow;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    sep[0]=';'; sep[1]=0;
    nsep = 1;
    quote = '"';
    escape = '\\';

    if (nrhs < 1 || nrhs > 4 || nlhs > 2) {
        mexErrMsgTxt("Usage: [numbers, text] = swallowcsv(filename, [sep, [quote, [escape]]])");
    }
    
    for (int i = 0; i < nrhs; i++) {
        if (!mxIsChar(prhs[i]) || mxIsEmpty(prhs[i])) {
            mexErrMsgTxt("Usage: [numbers, text] = swallowcsv(filename, [sep, [quote, [escape]]])");
        }
    }
    
    if (nrhs > 1) {
        nsep = mxGetNumberOfElements(prhs[1]);
        if (nsep<1 || nsep>80 || mxGetString(prhs[1], sep, nsep+1)) mexErrMsgTxt("sep is no valid string or too long");
        if (nrhs > 2) {
            quote = mxGetChars(prhs[2])[0] & 0xff;
            if (nrhs > 3) {
                escape = mxGetChars(prhs[3])[0] & 0xff;
            }
        }
    }
    
    char filename[1024];
    mxGetString(prhs[0], filename, 1024);
    FILE *f = fopen(filename, "rb");
    if (!f) {
        mexErrMsgTxt("Couldn't open specified file");
    }
    vector<vector<string> > data;
    size_t max_columns = 0;
    int nrows;
    
    parse_file(f, data, max_columns, nrows, NULL, NULL);
    
    mwSize dims[] = {(mwSize) nrows, (mwSize) max_columns};
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double *p = mxGetPr(plhs[0]);
    mxArray *pcell = NULL;
    if (nlhs==2) {plhs[1] = mxCreateCellArray(2, dims); pcell = plhs[1];}

    if (data.size()==nrows) {    
      for (int i = 0; i < data.size(); i++) {
          vector<string> &row(data[i]);
          for (int j = 0; j < row.size(); j++) 
            insert_elem(row[j].c_str(),row[j].size(),i,j,nrows,p,pcell);
          for (int j = row.size(); j < max_columns; j++)
            p[j * data.size() + i] = NaN;
      }
    } else {
      rewind(f);
      parse_file(f, data, max_columns, nrows, p, pcell);
    }
    fclose(f);      
}
