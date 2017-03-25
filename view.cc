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
    size_t chc = 0;
    for (int j = 1; j < argc; j++) {
        // Load GDML file with name
        G4GDMLParser p;
        G4cout << "Started reading (may take a while)..." << G4endl;
        p.Read(argv[j], false);
        G4cout << "Done reading..." << G4endl;
        GeoOption g;
        g.name = G4String(argv[j]);
        g.cons = NULL;
        g.cache = p.GetWorldVolume();
        opts.push_back(g);
        G4cout << "Done converting..." << G4endl;
        p.Clear();
    }

    Viewer v(opts, chc);
    return qapp.exec();
}
