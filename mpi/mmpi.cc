
#include <cassert>
#include <cstdio>
#include <iostream>
#include <unordered_map>
#include <string>

#include "Riostream.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "THnBase.h"
#include "TKey.h"
#include "TList.h"
#include "TMemFile.h"
#include "TTree.h"

using namespace std;

Int_t total_histos;
Int_t total_bins;
Int_t defined_bins;
Int_t total_histos_sparse;
Int_t total_bins_sparse;
Int_t defined_bins_sparse;

unordered_map<int, string> hist_names;
unordered_map<int, unordered_map<int, double> > hist;

static void Usage(const char *progname) {
  printf("usage: mpirun %s <ROOT histograms>\n", progname);
}

static string MkName(const TObject *o) {
  TString result = gDirectory->GetPath();
  result.Append("/");
  result.Append(o->GetName());
  return std::string(result.Data(), result.Length());
}

static void ProcessSparse(const THnBase *h) {
  total_histos_sparse++;
  total_bins_sparse += h->GetNbins();
}

static void ProcessHisto(const TH1 *h) {
  total_histos++;
  Int_t nx = h->GetNbinsX();
  Int_t ny = h->GetNbinsY();
  Int_t nz = h->GetNbinsZ();
  assert(nx > 0);
  assert(ny > 0);
  assert(nz > 0);
  defined_bins += nx*ny*nz;
  cerr << "    This histogram has bins " << nx << ", " << ny << ", " << nz << endl;

  string name = MkName(h);
  hist_names[total_histos] = name;
  hist[total_histos] = unordered_map<int, double>();

  Double_t all_sum = 0.0;
  for (Int_t x = 0; x <= nx; ++x) {
    if (ny > 1) {
      for (Int_t y = 0; y < ny; ++y) {
        if (nz > 1) {
          for (Int_t z = 0; z < nz; ++z) {
            Double_t this_bin = h->GetBinContent(x, y, z);
            if (this_bin > 0) {
              all_sum += this_bin;
              total_bins++;
              hist[total_histos][(nx+1)*(ny+1)*z + (nx+1)*y + x] = this_bin;
              //PPrintBin((nx+1)*(ny+1)*z + (nx+1)*y + x, this_bin);
            }
          }
        } else {
          Double_t this_bin = h->GetBinContent(x, y);
          if (this_bin > 0) {
            all_sum += this_bin;
            total_bins++;
            hist[total_histos][(nx+1)*y + x] = this_bin;
          }
        }
      }
    } else {
      Double_t this_bin = h->GetBinContent(x);
      if (this_bin > 0) {
        total_bins++;
        all_sum += this_bin;
        hist[total_histos][x] = this_bin;
      }
    }
  }
  cerr << "      All_sums " << all_sum << endl;
}

static void WalkRootfile(TFile *source) {
  TString path(gDirectory->GetPath());
  source->cd(path);
  cerr << "Browsing " << path << endl;

  // gain time, do not add the objects in the list in memory
  TH1::AddDirectory(kFALSE);

  // loop over all keys in this directory
  TIter nextkey(gDirectory->GetListOfKeys());
  TKey *key, *oldkey=0;
  while ((key = (TKey*)nextkey())) {
     //keep only the highest cycle number for each key
     if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

     // read object from first source file
     TObject *obj = key->ReadObj();

     cerr << "NAME IS " << obj->GetName() << endl;

     if (obj->IsA()->InheritsFrom(TH1::Class())) {
        cerr << "   HELLO, we have the histogram " << obj->GetName() << endl;
        TH1 *h = (TH1*)obj;
        ProcessHisto(h);
     } else if ( obj->IsA()->InheritsFrom( THnBase::Class())) {
        cerr << "   HELLO-SPARSE " << obj->GetName() << endl;
        THnBase *h = (THnBase*)obj;
        ProcessSparse(h);
     } else if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
        // it's a subdirectory
        //cerr << "Found subdirectory " << obj->GetName() << endl;
        source->cd( path + "/" + obj->GetName() );
        WalkRootfile(source);
        source->cd(path);
     } else if (obj->IsA()->InheritsFrom(TCollection::Class())) {
       TCollection *coll = (TCollection *)obj;
       //cerr << "List of something in " << obj->GetName() << endl;
       TIter nextelem(coll);
       TObject *elem;
       while ((elem = nextelem())) {
         if (elem->IsA()->InheritsFrom( TH1::Class() )) {
           cerr << "   HELLO, we have the histogram " << elem->GetName() << endl;
           TH1 *h = (TH1 *)elem;
           ProcessHisto(h);
         } else if (elem->IsA()->InheritsFrom( THnBase::Class() )) {
           cerr << "   HELLO-SPARSE " << elem->GetName() << endl;
           THnBase *h = (THnBase *)elem;
           ProcessSparse(h);
         }
       }
     } else {
        // object is of no type that we know or can handle
        cerr << "Unknown object type, name: "
          << obj->GetName() << " title: " << obj->GetTitle() << endl;
     }
     delete obj;
  }  // while (( TKey *key = (TKey*)nextkey()))
}

static void ReadHistos(char *input_path) {
  WalkRootfile(new TFile(input_path));
  cerr << "Total histos: " << total_histos << endl;
  cerr << "Total bins: " << total_bins << endl;
  cerr << "Defined bins: " << defined_bins << endl;
  cerr << "Total sparse histos: " << total_histos_sparse << endl;
  cerr << "Total sparse bins: " << total_bins_sparse << endl;
  cerr << "Defined sparse bins: " << defined_bins_sparse << endl;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    Usage(argv[0]);
    return 1;
  }
  ReadHistos(argv[1]);

  return 0;
}
