/*
 * The EDI Tools application.
 *
 * Copyright (C) 2024 Alexander Grayver <agrayver.geophysics@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "app/MainWindow.h"
#include <QApplication>
#include <QDir>
#include <QFileInfo>

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);

  QCoreApplication::setApplicationName("EDITools");

  // Portable distributions keep projection data and CA certificates beside
  // the executable (Windows) or one directory above it (Linux).
  const QDir executable(QCoreApplication::applicationDirPath());
  for(const auto &relative : {QStringLiteral("share"), QStringLiteral("../share")}) {
    const QDir data(executable.filePath(relative));
    if(!QFileInfo::exists(data.filePath("proj/proj.db"))) continue;
    qputenv("PROJ_DATA", data.filePath("proj").toUtf8());
    qputenv("GDAL_DATA", data.filePath("gdal").toUtf8());
    const auto certificates = data.filePath("certs/cacert.pem");
    if(QFileInfo::exists(certificates)) {
      qputenv("SSL_CERT_FILE", certificates.toUtf8());
      qputenv("CURL_CA_BUNDLE", certificates.toUtf8());
    }
    break;
  }

  MainWindow w;
  w.show();

  return a.exec();
}
