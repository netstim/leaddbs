
#ifndef STRAEMTHREAD_H
#define STRAEMTHREAD_H

#include <QMutex>
#include <QSize>
#include <QThread>
#include <QWaitCondition>


class qAlphaOmegaStatusThread : public QThread
{
  Q_OBJECT

public:

  qAlphaOmegaStatusThread(QObject *parent = nullptr);
  ~qAlphaOmegaStatusThread();

signals:

  void connectionStatusModified(bool* deviceIsConnected);
  void distanceToTargetModified(float *distanceToTargetMiliM);

protected:

  void run() override;

private:

  float GetDistanceToTargetMiliM();

  QMutex Mutex;
  QWaitCondition Condition;
  bool Abort = false;

};

#endif