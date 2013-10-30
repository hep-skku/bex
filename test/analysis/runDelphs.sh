#!/bin/bash

DELPHES=/users/jhgoh/work/Blackhole/Delphes/DelphesHepMC

$DELPHES delphesCMS.tcl delphes_ADD_10D_100TeV.root hadevent_ADD_10D_100TeV.hepmc &
$DELPHES delphesCMS.tcl delphes_ADD_10D_14TeV.root  hadevent_ADD_10D_14TeV.hepmc  &
$DELPHES delphesCMS.tcl delphes_ADD_10D_7TeV.root   hadevent_ADD_10D_7TeV.hepmc   &
$DELPHES delphesCMS.tcl delphes_ADD_10D_8TeV.root   hadevent_ADD_10D_8TeV.hepmc   &
$DELPHES delphesCMS.tcl delphes_RS_10D_100TeV.root  hadevent_RS_10D_100TeV.hepmc  &
$DELPHES delphesCMS.tcl delphes_RS_6D_100TeV.root   hadevent_RS_6D_100TeV.hepmc   &
$DELPHES delphesCMS.tcl delphes_RS_6D_14TeV.root    hadevent_RS_6D_14TeV.hepmc    &
$DELPHES delphesCMS.tcl delphes_RS_6D_8TeV.root     hadevent_RS_6D_8TeV.hepmc     &

wait

root -l -q -b 'delphes("delphes_ADD_10D_100TeV.root", "results_ADD_10D_100TeV.root")'
root -l -q -b 'delphes("delphes_ADD_10D_14TeV.root" , "results_ADD_10D_14TeV.root" )'
root -l -q -b 'delphes("delphes_ADD_10D_8TeV.root"  , "results_ADD_10D_8TeV.root"  )'
root -l -q -b 'delphes("delphes_ADD_10D_7TeV.root"  , "results_ADD_10D_7TeV.root"  )'
root -l -q -b 'delphes("delphes_RS_6D_100TeV.root"  , "results_RS_6D_100TeV.root"  )'
root -l -q -b 'delphes("delphes_RS_6D_14TeV.root"   , "results_RS_6D_14TeV.root"   )'
root -l -q -b 'delphes("delphes_RS_6D_8TeV.root"    , "results_RS_6D_8TeV.root"    )'

