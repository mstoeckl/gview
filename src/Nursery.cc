/* SPDX-License-Identifier: GPL-3.0-only */
#include "Nursery.hh"

#include <QThreadPool>

struct NData {
    QThreadPool pool;
};

Nursery::Nursery(int nthreads) {
    struct NData *d = new (struct NData)();
    data = d;
    d->pool.setMaxThreadCount(nthreads == 0 ? QThread::idealThreadCount()
                                            : nthreads);
}
Nursery::~Nursery() {
    struct NData *d = (struct NData *)data;
    d->pool.waitForDone();
    delete d;
}

void Nursery::startSoon(QRunnable *r, DeleteAction delrun) {
    struct NData *d = (struct NData *)data;
    switch (delrun) {
    case kUnset:
        break;
    case kDeleteAfter:
        r->setAutoDelete(true);
        break;
    case kKeepAfter:
        r->setAutoDelete(false);
        break;
    }
    d->pool.start(r);
}

int Nursery::idealThreadCount() const {
    struct NData *d = (struct NData *)data;
    return d->pool.maxThreadCount();
}
