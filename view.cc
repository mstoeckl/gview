#include "VectorWindow.hh"
#include "Viewer.hh"

#include <QApplication>
#include <QFile>
#include <QProcess>
#include <QTemporaryFile>
#include <QTextStream>

#include <G4GDMLParser.hh>

int usage() {
    fprintf(stderr,
            "Usage: view [--vector] [filename.gdml(.gz)]+ [trace.dat(.gz)]*\n");
    return -1;
}

int main(int argc, char **argv) {
    // Trick: we do basic init, and only use detector construction;
    // run an orthocamera (well, single point, and dx,dy,1 offset);
    // next up isa function handling point colors
    QApplication qapp(argc, argv);

    if (argc == 1) {
        return usage();
    }

    bool vector = false;
    std::vector<GeoOption> opts;
    std::vector<TrackData> tracks;
    for (int j = 1; j < argc; j++) {
        G4String fn(argv[j]);
        if (fn == "--vector") {
            vector = true;
            continue;
        }

        if (fn.size() >= 3 &&
            strcmp(fn.substr(fn.size() - 3, 3).data(), ".gz") == 0) {
            system("rm -f /tmp/copy.dat.gz");
            G4String ar = "cp " + fn + " /tmp/copy.dat.gz";
            system(ar.data());
            G4String gu = "gzip -df /tmp/copy.dat.gz";
            system(gu.data());
            fn = "/tmp/copy.dat";
        }

        if (fn.size() >= 5 &&
            strcmp(fn.substr(fn.size() - 5, 5).data(), ".gdml") == 0) {

            // Strip properties from file to avoid a GEANT bug
            G4cout << "Stripping extras from >" << fn << "<" << G4endl;
            QFile ifo(fn.data());
            QFile ofo("/tmp/cleaned.gdml");
            if (!ifo.open(QIODevice::ReadOnly) ||
                !ofo.open(QIODevice::WriteOnly)) {
                qFatal("Failed to open");
            }
            QTextStream in(&ifo);
            QTextStream out(&ofo);
            while (!in.atEnd()) {
                QString line = in.readLine();
                if (!line.contains("<property")) {
                    out << line << "\n";
                }
            }
            ifo.close();
            ofo.close();

            G4GDMLParser p;
            p.SetAddPointerToName(true);
            G4cout << "Started reading (may take a while)..." << G4endl;
            p.Read(ofo.fileName().toUtf8().constData(), false);
            G4cout << "Done reading..." << G4endl;
            GeoOption g;
            g.name = G4String(argv[j]);
            // Need to modify volume name to prevent collisions in lookup
            g.vol = p.GetWorldVolume();
            char buf[30];
            sprintf(buf, "-%d", j);
            G4String name = g.vol->GetName() + buf;
            g.vol->SetName(name);
            g.vol->GetLogicalVolume()->SetName(name);
            opts.push_back(g);
            G4cout << "Done converting..." << G4endl;
            p.Clear();
        } else if (fn.size() >= 4 &&
                   strcmp(fn.substr(fn.size() - 4, 4).data(), ".dat") == 0) {
            tracks.push_back(TrackData(fn.data()));
        } else {
            return usage();
        }
    }
    if (opts.size() == 0) {
        return usage();
    }

    if (!vector) {
        Viewer v(opts, tracks);
        return qapp.exec();
    } else {
        if (opts.size() != 1 || tracks.size() > 0) {
            return usage();
        }
        VectorWindow t(opts[0].name.c_str(), opts[0].vol);
        return qapp.exec();
    }
}
