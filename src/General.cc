#include "General.hh"

#include <QString>

QString FColor::hexName() const {
    QRgb color = rgba();
    int cr = qRed(color);
    int cg = qGreen(color);
    int cb = qBlue(color);
    const char *numbers = "0123456789abcdef";
    return QString("#%1%2%3%4%5%6")
        .arg(numbers[cr / 16])
        .arg(numbers[cr % 16])
        .arg(numbers[cg / 16])
        .arg(numbers[cg % 16])
        .arg(numbers[cb / 16])
        .arg(numbers[cb % 16]);
}
