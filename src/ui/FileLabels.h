#ifndef EDITOOLS_FILE_LABELS_H
#define EDITOOLS_FILE_LABELS_H

#include <QString>

namespace FileLabels {
// Display only the filename, including paths from projects saved on another OS.
inline QString fileName(QString path)
{
  path.replace('\\', '/');
  path = path.section('/', -1);
  if(path.size() >= 2 && path[0].isLetter() && path[1] == ':') path.remove(0, 2);
  return path;
}
}
#endif
