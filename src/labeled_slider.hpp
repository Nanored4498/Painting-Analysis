#include <QWidget>
#include <QSlider>
#include <QLabel>

class LabeledSlider: public QWidget {
Q_OBJECT

public:
	LabeledSlider(const QString &title, int min, int max, int val0=-1, double fact=1.0);
	~LabeledSlider();

signals:
	void valueChanged(int);

public slots:
	void changeValue(int);

private:
	QString title;
	QSlider *slider;
	QLabel *label;

	double fact=1.0;
};