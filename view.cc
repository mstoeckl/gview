#include <QApplication>

#include "Viewer.hh"
#include <G4GDMLParser.hh>

int usage() {
    fprintf(stderr, "Usage: view [filename.gdml]+\n");
    exit(1);
    return 1;
}

int main(int argc, char **argv) {
    // Trick: we do basic init, and only use detector construction;
    // run an orthocamera (well, single point, and dx,dy,1 offset);
    // next up isa function handling point colors
    QApplication qapp(argc, argv);

    if (argc == 1) {
        usage();
    }

    std::vector<GeoOption> opts;
    size_t chc = argc - 2;
    for (int j = 1; j < argc; j++) {
        // Load GDML file with name
        G4GDMLParser p;
        p.SetAddPointerToName(true);
        G4cout << "Started reading (may take a while)..." << G4endl;
        p.Read(argv[j], false);
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
    }

    Viewer v(opts, chc);
    return qapp.exec();
}
