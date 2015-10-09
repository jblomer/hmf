
#include <mpi.h>
#include <pthread.h>

#include <cassert>
#include <cstdio>
#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>

#include "murmur.h"

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

int world_rank, world_size;

Int_t total_histos;
Int_t total_bins;
Int_t defined_bins;
Int_t total_histos_sparse;
Int_t total_bins_sparse;
Int_t defined_bins_sparse;

unordered_map<int, string> hist_names;
unordered_map<int, unordered_map<int, double> > hist;

struct Bin {
  Bin() : nhist(0), nbin(0), val(.0), err(.0) { }
  Bin(int n, int b, double v) : nhist(n), nbin(b), val(v), err(.0) { }
  int nhist;
  int nbin;
  double val;
  double err;
};

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

static void *main_recv(void *data) {
  for (int i = 0; i < world_size; ++i) {
    unsigned size;
    MPI_Recv(&size, sizeof(size), MPI_BYTE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    Bin *recv_buf = (Bin *)malloc(size * sizeof(Bin));
    MPI_Recv(recv_buf, size * sizeof(Bin), MPI_BYTE, MPI_ANY_SOURCE, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("  I'm %d, got %dth data packet\n", world_rank, i);
    for (unsigned i = 0; i < size; ++i) {
      hist[recv_buf[i].nhist][recv_buf[i].nbin] += recv_buf[i].val;
    }
    free(recv_buf);
  }
  return NULL;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    Usage(argv[0]);
    return 1;
  }

  int claimed, provided;
  MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &provided);
  MPI_Query_thread(&claimed);
  printf("Query thread level= %d  Init_thread level= %d\n", claimed, provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  ReadHistos(argv[1]);

  MPI_Barrier(MPI_COMM_WORLD);

  printf("rank %d, filling dispatch vector\n", world_rank);

  vector<vector<Bin>> dispatch;
  for (int i = 0; i < world_size; ++i) {
    dispatch.push_back(vector<Bin>());
  }

  for (auto i = hist.begin(); i != hist.end(); ++i) {
    for (auto j = i->second.begin(); j != i->second.end(); ++j) {
      struct {
        int a;
        int b;
      } global_bin;
      global_bin.a = i->first;
      global_bin.b = j->second;
      uint32_t hashed = MurmurHash2(&global_bin, sizeof(global_bin), 42);
      dispatch[hashed % world_rank].push_back(
        Bin(i->first, j->first, j->second));
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  printf("rank %d, sending and receiving data\n", world_rank);

  pthread_t thread_recv;
  pthread_create(&thread_recv, NULL, main_recv, NULL);
  for (int i = 0; i < world_size; ++i) {
    printf("I'm rank %d, sending data to rank %d\n", world_rank, i);
    unsigned size = dispatch[i].size();
    MPI_Ssend(&size, sizeof(unsigned), MPI_BYTE, i, 0, MPI_COMM_WORLD);
    MPI_Ssend(&dispatch[i][0], dispatch[i].size() * sizeof(Bin), MPI_BYTE,
              i, 0, MPI_COMM_WORLD);
  }
  pthread_join(thread_recv, NULL);

  MPI_Barrier(MPI_COMM_WORLD);

  if (world_rank != 0) {
    MPI_Finalize();
    return 0;
  }

  printf("Here's your master thread, reading out results.\n");
  MPI_Finalize();
  return 0;
}
