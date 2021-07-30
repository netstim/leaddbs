#include "qAlphaOmegaStatusThread.h"

// STD Includes
#include <cmath> // NAN
// #include <stdlib.h> 

// Windows
#include <Windows.h>

// AlphaOmega SDK
#include "AOSystemAPI.h"
#include "AOTypes.h"
#include "StreamFormat.h"


static const float DRIVE_ZERO_POSITION_MILIM = 25.0;
static const int SLEEP_TIME_MILIS = 100;


qAlphaOmegaStatusThread::qAlphaOmegaStatusThread(QObject *parent)
    : QThread(parent)
{

}

qAlphaOmegaStatusThread::~qAlphaOmegaStatusThread()
{
  this->Mutex.lock();
  this->Abort = true;
  this->Condition.wakeOne();
  this->Mutex.unlock();
  wait();
}


float qAlphaOmegaStatusThread::GetDistanceToTargetMiliM()
{
	int32 nDepthUm = 0;
	EAOResult eAORes = (EAOResult)GetDriveDepth(&nDepthUm);

	if (eAORes == eAO_OK)
	{
    return DRIVE_ZERO_POSITION_MILIM - nDepthUm / 1000.0;
  }
  else
  {
    return NAN;
  }
}


void qAlphaOmegaStatusThread::run()
{

  bool deviceWasConnected = false;
  bool deviceIsConnected = false;

  float previousDistanceToTargetMiliM = NAN;
  float distanceToTargetMiliM = NAN;


  forever
  {

    Sleep(SLEEP_TIME_MILIS);

    this->Mutex.lock();

    deviceIsConnected = (isConnected() == eAO_CONNECTED);
    if (deviceIsConnected != deviceWasConnected)
    {
      emit connectionStatusModified(&deviceIsConnected);
      deviceWasConnected = deviceIsConnected;
    }


    distanceToTargetMiliM = this->GetDistanceToTargetMiliM();
    // distanceToTargetMiliM = 10.0 * ((float) rand()) / (float) RAND_MAX;
    if (distanceToTargetMiliM != previousDistanceToTargetMiliM)
    {
      emit distanceToTargetModified(&distanceToTargetMiliM);
      previousDistanceToTargetMiliM = distanceToTargetMiliM;
    }

    this->Mutex.unlock();

  }
}
