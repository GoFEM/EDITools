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

#ifndef MT_DATA_PLOT_H
#define MT_DATA_PLOT_H

#include "qcustomplot.h"

#include "MTStationData.h"

class MyCustomPlot: public QCustomPlot
{
public:
  MyCustomPlot(QWidget *parent):
    QCustomPlot(parent)
  {}

  void  keyPressEvent(QKeyEvent *event)
  {
    if(event->key() == Qt::Key_X)
      axisRect()->setRangeZoomAxes(nullptr, this->yAxis);
    else if(event->key() == Qt::Key_D)
      setSelectionRectMode(QCP::srmNone);

    QCustomPlot::keyPressEvent(event);
  }

  void keyReleaseEvent(QKeyEvent *event)
  {
    if(event->key() == Qt::Key_X)
      axisRect()->setRangeZoomAxes(this->xAxis, this->yAxis);
    else if(event->key() == Qt::Key_D)
      setSelectionRectMode(QCP::srmSelect);

    QCustomPlot::keyReleaseEvent(event);
  }
};

class MTDataPlot: public QObject
{
  Q_OBJECT

public:
  MTDataPlot(QCustomPlot* plot);

  virtual void set_observed_data(MTStationData &data, bool rescaleAxes = true) = 0;
  virtual void set_predicted_data(const MTStationData &data, bool rescaleAxes = false) = 0;

  void set_associated_plot(MTDataPlot &plot);

  void set_error_bars_visible(bool on);
  void set_masking_mode(bool on);

  std::string get_graph_data_type_name(const QCPGraph *graph) const;

public slots:
  void maskSelectedData();
  void invMaskSelectedData();
  void unmaskSelectedData();
  void dataSelected(bool selected);
  void plotContextRequest(QPoint pos);
  void showPointToolTip(QMouseEvent* event);

protected:
  void set_graph_data(const std::vector<std::vector<bool>> &mask,
                      const std::vector<std::vector<double>> &data,
                      const std::vector<std::vector<double>> &data_err,
                      const std::vector<double> &frequencies);

  void set_graph_responses(const std::vector<std::vector<double>> &data,
                           const std::vector<double> &frequencies);

  void set_layout_generic(const std::vector<QString> &data_graph_names,
                          const std::vector<QString> &masked_graph_names);

  void mask_selected_data(bool on);

protected:
  MTDataPlot* m_associated_plot;
  QCustomPlot* m_plot;
  QMenu* m_contextMenu;

  MTStationData *m_data;

  std::vector<QCPErrorBars*> m_errorBars;

  std::map<std::string, RealDataType> m_name2type;
};

#endif
