#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "include/RhoPlusOverlay.h"
#include "include/StationListDelegate.h"
#include "include/TFViewControls.h"
#include <QSignalBlocker>
#include <algorithm>
#include <cmath>

TFViewRange MainWindow::viewRange(int index) const
{
  if(index == 0) {
    const auto range = plot11->xAxis->range();
    return {periodRange.automatic, range.lower, range.upper};
  }
  const int plotIndex = (index - 1) % 2;
  const bool diagonal = index >= 3;
  const auto range = (plotIndex == 0 ? plot11 : plot12)->axisRect(diagonal ? 1 : 0)->axis(QCPAxis::atLeft)->range();
  return {plotHandlers[plotIndex]->y_axis_autoscale(diagonal), range.lower, range.upper};
}

void MainWindow::setViewRange(int index, TFViewRange range)
{
  if(index == 0) {
    periodRange = range;
    if(range.automatic) fitViewRange(index);
    else plot11->xAxis->setRange(range.lower, range.upper);
  } else {
    const int plotIndex = (index - 1) % 2;
    const bool diagonal = index >= 3;
    plotHandlers[plotIndex]->set_y_axis_autoscale(range.automatic, diagonal);
    plotHandlers[plotIndex]->set_y_axis_range(range.lower, range.upper, diagonal);
    if(range.automatic) fitViewRange(index);
    else (plotIndex == 0 ? plot11 : plot12)->axisRect(diagonal ? 1 : 0)->axis(QCPAxis::atLeft)->setRange(range.lower, range.upper);
  }
  for(auto *plot: {plot11, plot12, plot21, plot22}) plot->replot(QCustomPlot::rpQueuedReplot);
}

void MainWindow::fitViewRange(int index)
{
  if(index == 0) {
    const auto *item = ui->stationList->currentItem();
    if(!mtSurvey || !item) return;
    const auto &frequencies = mtSurvey->get_station_data(item->text().toStdString()).frequencies();
    if(frequencies.empty()) return;
    const auto bounds = std::minmax_element(frequencies.begin(), frequencies.end());
    if(*bounds.first > 0. && std::isfinite(*bounds.second))
      plot11->xAxis->setRange(.5 / *bounds.second, 2. / *bounds.first);
  } else {
    const bool phase = index == 2 || index == 4;
    auto *axis = (phase ? plot12 : plot11)->axisRect(index >= 3 ? 1 : 0)->axis(QCPAxis::atLeft);
    if(phase) {
      if(index == 2 && ui->actionPhase_wrap->isChecked()) axis->setRange(0, 90);
      else axis->setRange(-180, 180);
    } else axis->rescale(true);
  }
}

void MainWindow::applyViewStyle()
{
  for(auto &plot: plotHandlers) plot->set_view_style(viewControls->style());
  rhoPlusOverlay->set_colors(viewControls->style().colors);
}

void MainWindow::updateStationDetails()
{
  if(!mtSurvey) return;
  const QSignalBlocker blocker(ui->stationList);
  for(int i = 0; i < ui->stationList->count(); ++i) {
    auto *item = ui->stationList->item(i);
    const auto count = mtSurvey->get_station_data(item->text().toStdString()).frequencies().size();
    item->setData(StationListDelegate::PeriodCount, qulonglong(count));
    const auto details = tr("%1\n%2 periods\n%3\n%4").arg(item->text()).arg(count)
      .arg(item->data(StationListDelegate::Coordinates).toString(), item->data(StationListDelegate::CoordinateReference).toString());
    item->setToolTip(details);
    item->setData(Qt::AccessibleTextRole, details);
  }
}
