#include "Pythia.h"
#include "HepMCInterface.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/Units.h"

#include <string>

using namespace std;
using namespace Pythia8;

int usage()
{
  cout << "runPythia - Read LHE file and launch pythia8 to do hadronization" << endl;
  cout << "Usage : runPythia INPUTFILE OUTPUTFILE" << endl;
  return 1;
}

int main(int argc, char* argv[])
{
  if ( argc < 3 ) return usage();

  string inputFileName = argv[1];
  string outputFileName = argv[2];

  Pythia pythia;
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = " + inputFileName);
  pythia.readString("PartonLevel:all = true");
  pythia.init();

  HepMC::I_Pythia8 hepMC;
  HepMC::IO_GenEvent ascii_io(outputFileName, std::ios::out);

  for ( int i=0; ; ++i )
  {
    if ( !pythia.next() and pythia.info.atEndOfFile() ) break;

    HepMC::GenEvent* genEvent = new HepMC::GenEvent();
    hepMC.fill_next_event(pythia, genEvent);

    ascii_io << genEvent;
    delete genEvent;
  }

  pythia.stat();

  return 0;
}

