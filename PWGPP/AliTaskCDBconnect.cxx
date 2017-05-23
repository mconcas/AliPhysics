/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliTaskCDBconnect.h"

#include <TChain.h>
#include <TFile.h>
#include <TRegexp.h>
#include <TGeoGlobalMagField.h>
#include "TGeoManager.h"
 
#include "AliAnalysisManager.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliVEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"

ClassImp(AliTaskCDBconnect)

//______________________________________________________________________________
AliTaskCDBconnect::AliTaskCDBconnect():
  AliAnalysisTask(),
  fFallBackToRaw(kFALSE),
  fRun(0),
  fLock(0),
  fStorage(),
  fSpecCDBUri(),
  fGRPManager(NULL)
{
  // Dummy constructor
  fSpecCDBUri.SetOwner();
}

//______________________________________________________________________________
AliTaskCDBconnect::AliTaskCDBconnect(const char* name, const char *storage, Int_t run, Bool_t fallback)
  :AliAnalysisTask(name, "ESD analysis tender car"),
   fFallBackToRaw(fallback),
   fRun(run),
   fLock(0),
   fStorage(storage),
   fSpecCDBUri(),
   fGRPManager(NULL)
{
  // Default constructor
  fSpecCDBUri.SetOwner();
  DefineInput (0, TChain::Class());
  if (run>0) InitGRP();
}

//______________________________________________________________________________
AliTaskCDBconnect::~AliTaskCDBconnect()
{
  // Destructor
  delete fGRPManager;
}  

//______________________________________________________________________________
void AliTaskCDBconnect::InitGRP()
{
  // Initialize geometry and mag. field
  AliCDBManager *cdb = AliCDBManager::Instance();
  Bool_t useCVMFS = kFALSE;
  TString storage = fStorage.Strip(TString::kTrailing,'/');
  TObjArray *inpStor = storage.Tokenize("[:]");

  // If CVMFS -> set storage at chosen path
  if (((TObjString *)(inpStor->At(0)))->String() == "cvmfs") {
    useCVMFS = kTRUE;
    if (inpStor->GetEntries()<3 && inpStor->GetEntries()>1) {
      gSystem->Setenv("OCDB_PATH", ((TObjString *)(inpStor->At(1)))->String());
    } 
    else if (inpStor->GetEntries() == 1) {
      gSystem->Setenv("OCDB_PATH", "/cvmfs/alice-ocdb.cern.ch/");
    }
    AliInfoF("Setting default storage to %s",fStorage.Data());
    cdb->SetDefaultStorage("raw://"); // "raw" to "cvmfs" switch is transparently managed by AliCDBManager
  } else {
    AliInfoF("Setting default storage to %s",fStorage.Data());
    cdb->SetDefaultStorage(fStorage);
  }

  // set specific storages
  for (Int_t i = 0; i < fSpecCDBUri.GetEntriesFast(); i++) {
    TNamed* obj = (TNamed*)fSpecCDBUri[i];
    if (!obj) continue;
    UInt_t vsv = obj->GetUniqueID();
    Int_t ver    = int(vsv>>16)-1;
    Int_t subver = int(vsv&0xffff)-1;
    cdb->SetSpecificStorage(obj->GetName(), obj->GetTitle(), ver, subver);
  }
  if (cdb->GetRun()!=fRun) {    
    fLock = cdb->SetLock(kFALSE,fLock);
    cdb->SetRun(fRun);
    fLock = cdb->SetLock(kTRUE,fLock);
  }
  if (gSystem->AccessPathName("OCDB.root",kFileExists)==0) cdb->SetSnapshotMode("OCDB.root"); 
  //
  if (!fGRPManager) fGRPManager = new AliGRPManager();
  AliInfo("AliCDBconnect: #### Loading GRP to init B-field...");
  if(!fGRPManager->ReadGRPEntry()) AliFatal("Cannot get GRP entry"); 
  if(!fGRPManager->SetMagField())  AliFatal("Problem with magnetic field setup"); 
  //
  // geometry
  if (!gGeoManager) {
    AliInfo("AliCDBconnect: #### Loading geometry...");
    AliGeomManager::LoadGeometry("geometry.root");
    if(!AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC TRD TOF EMCAL PHOS MUON")) AliWarning("Problem with align objects");
  }  
}

//______________________________________________________________________________
void AliTaskCDBconnect::CreateOutputObjects()
{
  // Init CDB locally if run number is defined.
  //
  //  try to init before the analysis set
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) AliFatal("No analysis manager");
  if (fRun>0 && !fGRPManager) {
    // in the proof or plugin mode the initialization done in the constructor is not available
    InitGRP();
  }
  else {
    AliInfo("Run number is not available at this stage, InitGRP will be called in the execution loop");
  }
}

//______________________________________________________________________________
void AliTaskCDBconnect::ConnectInputData(Option_t* option)
{
  // Connect the input data, create CDB manager.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) AliFatal("No analysis manager");  
  AliAnalysisTask::ConnectInputData(option);
  Int_t run = AliAnalysisManager::GetAnalysisManager()->GetRunFromPath();
  if (run<=0) {
    AliWarning("AliTaskCDBconnect: Could not set run from path");
    return;
  }
  if (fRun != run) {
    fRun = run;
    InitGRP();
  }
}


//______________________________________________________________________________
void AliTaskCDBconnect::Exec(Option_t* /*option*/)
{
  // Execute all supplied analysis of one event. Notify run change via RunChanged().
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) AliFatal("No analysis manager");
  AliInputEventHandler* inp = (AliInputEventHandler*)mgr->GetInputEventHandler();
  if (!inp) AliFatal("No input event handler connected");
  AliVEvent* ev = inp->GetEvent();
  if (!ev) AliFatal("No event returned");
  int run = ev->GetRunNumber();
  // Intercept when the run number changed
  if (fRun != run) {
    fRun = run;
    InitGRP();
  }
}

//______________________________________________________________________________
void AliTaskCDBconnect::SetSpecificStorage(const char* calibType, const char* dbString, Int_t version, Int_t subVersion)
{
  // Set a specific storage
  TNamed *nmpath = new TNamed(calibType,dbString);
  if (version<0) version = -1;
  if (subVersion<0) subVersion = -1;
  nmpath->SetUniqueID((UInt_t(version+1)<<16)+UInt_t(subVersion+1));
  fSpecCDBUri.AddLast(nmpath);
}
