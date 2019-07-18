#include "labeled_slider.hpp"

#include <QVBoxLayout>

LabeledSlider::LabeledSlider(const QString &title, int min, int max, int val0, double p_fact):
	title(title), fact(p_fact), QWidget() {
	
	QVBoxLayout *layout = new QVBoxLayout();
	layout->setContentsMargins(0, 0, 0, 0);
	layout->setSpacing(0);
	setLayout(layout);

	if(val0 == -1) val0 = min;
	max = qMax(min, max);
	val0 = qMin(max, qMax(min, val0));

	label = new QLabel(title + " (" + QString::number(val0*fact) + ")");
	layout->addWidget(label);
	
	slider = new QSlider(Qt::Orientation::Horizontal);
	slider->setRange(min, max);
	slider->setValue(val0);
	layout->addWidget(slider);
	connect(slider, SIGNAL(valueChanged(int)), this, SLOT(changeValue(int)));
}

LabeledSlider::~LabeledSlider() {
	delete label;
	delete slider;
}

void LabeledSlider::changeValue(int v) {
	label->setText(title + " (" + QString::number(v*fact) + ")");
	emit valueChanged(v);
}