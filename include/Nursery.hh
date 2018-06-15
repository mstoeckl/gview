/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

class QRunnable;

typedef enum { kUnset, kDeleteAfter, kKeepAfter } DeleteAction;

/**
 * Task nursery. Bbecause parentless threads of any type are evil.
 * Might be called QScopedThreadPool.
 */
class Nursery {
public:
    Nursery(int nthreads = 0);
    ~Nursery();

    void startSoon(QRunnable *r, DeleteAction delrun = kUnset);
    int idealThreadCount() const;

private:
    void *data;
};

/**
 * To signify when to stop.
 */
class CancellationToken {
public:
    CancellationToken() : canceled(false) {}
    volatile bool canceled;

private:
    CancellationToken(const CancellationToken &) = delete;
    CancellationToken &operator=(const CancellationToken &) = delete;
};
