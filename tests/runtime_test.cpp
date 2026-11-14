#include <QApplication>
#include <QOperatingSystemVersion>
#include <QSslSocket>
#include <QStyle>
#include <QStyleFactory>
#include <iostream>

int main(int argc, char **argv)
{
#ifdef Q_OS_WIN
  // Exercise the native style even when the other GUI tests run offscreen.
  qputenv("QT_QPA_PLATFORM", "windows");
  qunsetenv("QT_STYLE_OVERRIDE");
#endif
  QApplication app(argc, argv);
  std::cout << "Qt " << qVersion() << "; style "
            << app.style()->objectName().toStdString() << '\n';
#ifdef Q_OS_WIN
  if(!QStyleFactory::keys().contains("windows11", Qt::CaseInsensitive)) {
    std::cerr << "Native Windows style plugin is missing\n";
    return 1;
  }
  if(QOperatingSystemVersion::current() >= QOperatingSystemVersion::Windows11 &&
     app.style()->objectName().compare("windows11", Qt::CaseInsensitive) != 0) {
    std::cerr << "Windows 11 did not select its native style\n";
    return 1;
  }
#endif
  if(!QSslSocket::supportsSsl()) {
    std::cerr << "TLS backend is unavailable\n";
    return 1;
  }
  return 0;
}
