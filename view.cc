#include <QApplication>

#include "Viewer.hh"
#include <G4GDMLParser.hh>

int usage() {
    fprintf(stderr, "Usage: view [filename.gdml]+ [trace.dat]*\n");
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

    std::vector<GeoOption> opts;
    std::vector<TrackData> tracks;
    for (int j = 1; j < argc; j++) {
        G4String fn(argv[j]);
        if (fn.size() >= 3 &&
            strcmp(fn.substr(fn.size() - 3, 3).data(), ".gz") == 0) {
            // Undo gzip!
            G4String ar = "cp " + fn + " /tmp/copy.dat.gz";
            system(ar.data());
            G4String gu = "gzip -df /tmp/copy.dat.gz";
            system(gu.data());
            fn = "/tmp/copy.dat";
        }

        if (fn.size() >= 5 &&
            strcmp(fn.substr(fn.size() - 5, 5).data(), ".gdml") == 0) {
            G4GDMLParser p;
            p.SetAddPointerToName(true);
            G4cout << "Started reading (may take a while)..." << G4endl;
            p.Read(fn, false);
            G4cout << "Done reading..." << G4endl;
            GeoOption g;
            g.name = G4String(argv[j]);
            g.cons = NULL;
            // Need to modify volume name to prevent collisions in lookup
            g.cache = p.GetWorldVolume();
            char buf[30];
            sprintf(buf, "-%d", j);
            G4String name = g.cache->GetName() + buf;
            g.cache->SetName(name);
            g.cache->GetLogicalVolume()->SetName(name);
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

    Viewer v(opts, tracks);
    return qapp.exec();
}
