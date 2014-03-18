/*
 *  This file is part of RawTherapee.
 */
#include "gradient.h"
#include "rtimage.h"
#include <iomanip>
#include "../rtengine/rt_math.h"

using namespace rtengine;
using namespace rtengine::procparams;

Gradient::Gradient () : Gtk::VBox(), FoldableToolPanel(this), EditSubscriber(ET_OBJECTS), lastObject(-1), draggedPointOldAngle(-1000.)
{
	set_border_width(4);

	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabled->set_active (false);
	enaConn  = enabled->signal_toggled().connect( sigc::mem_fun(*this, &Gradient::enabledChanged) );

	edit = Gtk::manage (new Gtk::ToggleButton());
	edit->add (*Gtk::manage (new RTImage ("editmodehand.png")));
	edit->set_tooltip_text(M("EDIT_OBJECT_TOOLTIP"));
	editConn = edit->signal_toggled().connect( sigc::mem_fun(*this, &Gradient::editToggled) );

	strength = Gtk::manage (new Adjuster (M("TP_GRADIENT_STRENGTH"), -5, 5, 0.01, 0));
	strength->set_tooltip_text (M("TP_GRADIENT_STRENGTH_TOOLTIP"));
	strength->setAdjusterListener (this);

	degree = Gtk::manage (new Adjuster (M("TP_GRADIENT_DEGREE"), -180, 180, 1, 0));
	degree->set_tooltip_text (M("TP_GRADIENT_DEGREE_TOOLTIP"));
	degree->setAdjusterListener (this);

	feather = Gtk::manage (new Adjuster (M("TP_GRADIENT_FEATHER"), 0, 100, 1, 25));
	feather->set_tooltip_text (M("TP_GRADIENT_FEATHER_TOOLTIP"));
	feather->setAdjusterListener (this);

	centerX = Gtk::manage (new Adjuster (M("TP_GRADIENT_CENTER_X"), -100, 100, 1, 0));
	centerX->set_tooltip_text (M("TP_GRADIENT_CENTER_X_TOOLTIP"));
	centerX->setAdjusterListener (this);

	centerY = Gtk::manage (new Adjuster (M("TP_GRADIENT_CENTER_Y"), -100, 100, 1, 0));
	centerY->set_tooltip_text (M("TP_GRADIENT_CENTER_Y_TOOLTIP"));
	centerY->setAdjusterListener (this);

	Gtk::HBox* enaBox = Gtk::manage (new Gtk::HBox());
	enaBox->pack_start(*enabled);
	enaBox->pack_end(*edit, false, false, 0);
	pack_start(*enaBox);
	pack_start(*Gtk::manage (new  Gtk::HSeparator()), Gtk::PACK_EXPAND_WIDGET, 4);
	pack_start (*strength);
	pack_start (*degree);
	pack_start (*feather);
	pack_start (*centerX);
	pack_start (*centerY);

	// Instantiating the Editing geometry; positions will be initialized later
	Line *hLine, *vLine, *featherLine[2];
	Circle *centerCircle;

	// Visible geometry
	hLine = new Line(); hLine->innerLineWidth=2;
	vLine = new Line();
	hLine->datum = vLine->datum = Geometry::IMAGE;

	featherLine[0] = new Line();  featherLine[0]->innerLineWidth=2;
	featherLine[1] = new Line();  featherLine[1]->innerLineWidth=2;
	featherLine[0]->datum = featherLine[1]->datum = Geometry::IMAGE;

	centerCircle = new Circle();
	centerCircle->datum = Geometry::IMAGE;
	centerCircle->radiusInImageSpace = false;
	centerCircle->radius = 6;
	centerCircle->filled = true;

	EditSubscriber::visibleGeometry.push_back( hLine );
	EditSubscriber::visibleGeometry.push_back( vLine );
	EditSubscriber::visibleGeometry.push_back( featherLine[0] );
	EditSubscriber::visibleGeometry.push_back( featherLine[1] );
	EditSubscriber::visibleGeometry.push_back( centerCircle );

	// MouseOver geometry
	hLine = new Line(); hLine->innerLineWidth=2;
	vLine = new Line();
	hLine->datum = vLine->datum = Geometry::IMAGE;

	featherLine[0] = new Line();  featherLine[0]->innerLineWidth=2;
	featherLine[1] = new Line();  featherLine[1]->innerLineWidth=2;
	featherLine[0]->datum = featherLine[1]->datum = Geometry::IMAGE;

	centerCircle = new Circle();
	centerCircle->datum = Geometry::IMAGE;
	centerCircle->radiusInImageSpace = false;
	centerCircle->radius = 30;
	centerCircle->filled = true;

	EditSubscriber::mouseOverGeometry.push_back( hLine );
	EditSubscriber::mouseOverGeometry.push_back( vLine );
	EditSubscriber::mouseOverGeometry.push_back( featherLine[0] );
	EditSubscriber::mouseOverGeometry.push_back( featherLine[1] );
	EditSubscriber::mouseOverGeometry.push_back( centerCircle );

	show_all();
}

Gradient::~Gradient() {
	for (std::vector<Geometry*>::const_iterator i = visibleGeometry.begin(); i != visibleGeometry.end(); ++i) {
		delete *i;
	}
	for (std::vector<Geometry*>::const_iterator i = mouseOverGeometry.begin(); i != mouseOverGeometry.end(); ++i) {
		delete *i;
	}
}

void Gradient::read (const ProcParams* pp, const ParamsEdited* pedited)
{
	disableListener ();

	if (pedited) {
		degree->setEditedState (pedited->gradient.degree ? Edited : UnEdited);
		feather->setEditedState (pedited->gradient.feather ? Edited : UnEdited);
		strength->setEditedState (pedited->gradient.strength ? Edited : UnEdited);
		centerX->setEditedState (pedited->gradient.centerX ? Edited : UnEdited);
		centerY->setEditedState (pedited->gradient.centerY ? Edited : UnEdited);
		enabled->set_inconsistent (multiImage && !pedited->gradient.enabled);
	}

	enaConn.block (true);
	enabled->set_active (pp->gradient.enabled);
	enaConn.block (false);
	degree->setValue (pp->gradient.degree);
	feather->setValue (pp->gradient.feather);
	strength->setValue (pp->gradient.strength);
	centerX->setValue (pp->gradient.centerX);
	centerY->setValue (pp->gradient.centerY);

	lastEnabled = pp->gradient.enabled;

	updateGeometry (pp->gradient.centerX, pp->gradient.centerY, pp->gradient.feather, pp->gradient.degree);

	enableListener ();
}

void Gradient::updateGeometry(int centerX_, int centerY_, double feather_, double degree_) {
	EditDataProvider* dataProvider = getEditProvider();
	if (dataProvider) {
		int imW, imH;
		PolarCoord polCoord1, polCoord2;
		dataProvider->getImageSize(imW, imH);
		double decay = feather_ * sqrt(double(imW)*double(imW)+double(imH)*double(imH)) / 200.;

		Coord origin(imW/2+centerX_*imW/200.f, imH/2+centerY_*imH/200.f);

		Line *currLine;
		Circle *currCircle;
		// update horizontal line
		currLine = static_cast<Line*>(visibleGeometry.at(0));
		polCoord1.set(1500.f, float(-degree_+180));  currLine->begin.setFromPolar(polCoord1);  currLine->begin += origin;
		polCoord1.set(1500.f, float(-degree_    ));  currLine->end.setFromPolar  (polCoord1);  currLine->end   += origin;
		currLine = static_cast<Line*>(mouseOverGeometry.at(0));
		polCoord1.set(1500.f, float(-degree_+180));  currLine->begin.setFromPolar(polCoord1);  currLine->begin += origin;
		polCoord1.set(1500.f, float(-degree_    ));  currLine->end.setFromPolar  (polCoord1);  currLine->end   += origin;
		// update vertical line
		currLine = static_cast<Line*>(visibleGeometry.at(1));
		polCoord1.set( 700.f, float(-degree_+90 ));  currLine->begin.setFromPolar(polCoord1);  currLine->begin += origin;
		polCoord1.set( 700.f, float(-degree_+270));  currLine->end.setFromPolar  (polCoord1);  currLine->end   += origin;
		currLine = static_cast<Line*>(mouseOverGeometry.at(1));
		polCoord1.set( 700.f, float(-degree_+90 ));  currLine->begin.setFromPolar(polCoord1);  currLine->begin += origin;
		polCoord1.set( 700.f, float(-degree_+270));  currLine->end.setFromPolar  (polCoord1);  currLine->end   += origin;
		// update upper feather line
		currLine = static_cast<Line*>(visibleGeometry.at(2));
		polCoord2.set(decay, float(-degree_+270));
		polCoord1.set(350.f, float(-degree_+180));  currLine->begin.setFromPolar(polCoord1+polCoord2);  currLine->begin += origin;
		polCoord1.set(350.f, float(-degree_    ));  currLine->end.setFromPolar  (polCoord1+polCoord2);  currLine->end   += origin;
		currLine = static_cast<Line*>(mouseOverGeometry.at(2));
		polCoord1.set(350.f, float(-degree_+180));  currLine->begin.setFromPolar(polCoord1+polCoord2);  currLine->begin += origin;
		polCoord1.set(350.f, float(-degree_    ));  currLine->end.setFromPolar  (polCoord1+polCoord2);  currLine->end   += origin;
		// update lower feather line
		currLine = static_cast<Line*>(visibleGeometry.at(3));
		polCoord2.set(decay, float(-degree_+90));
		polCoord1.set(350.f, float(-degree_+180));  currLine->begin.setFromPolar(polCoord1+polCoord2);  currLine->begin += origin;
		polCoord1.set(350.f, float(-degree_    ));  currLine->end.setFromPolar  (polCoord1+polCoord2);  currLine->end   += origin;
		currLine = static_cast<Line*>(mouseOverGeometry.at(3));
		polCoord1.set(350.f, float(-degree_+180));  currLine->begin.setFromPolar(polCoord1+polCoord2);  currLine->begin += origin;
		polCoord1.set(350.f, float(-degree_    ));  currLine->end.setFromPolar  (polCoord1+polCoord2);  currLine->end   += origin;
		// update circle's position
		currCircle = static_cast<Circle*>(visibleGeometry.at(4));
		currCircle->center = origin;
		currCircle = static_cast<Circle*>(mouseOverGeometry.at(4));
		currCircle->center = origin;
	}
}

void Gradient::write (ProcParams* pp, ParamsEdited* pedited)
{
	pp->gradient.degree = degree->getValue ();
	pp->gradient.feather = feather->getIntValue ();
	pp->gradient.strength = strength->getValue ();
	pp->gradient.centerX = centerX->getIntValue ();
	pp->gradient.centerY = centerY->getIntValue ();
	pp->gradient.enabled = enabled->get_active();

	if (pedited) {
		pedited->gradient.degree = degree->getEditedState ();
		pedited->gradient.feather = feather->getEditedState ();
		pedited->gradient.strength = strength->getEditedState ();
		pedited->gradient.centerX = centerX->getEditedState ();
		pedited->gradient.centerY = centerY->getEditedState ();
		pedited->gradient.enabled = !enabled->get_inconsistent();
	}
}

void Gradient::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{
	degree->setDefault (defParams->gradient.degree);
	feather->setDefault (defParams->gradient.feather);
	strength->setDefault (defParams->gradient.strength);
	centerX->setDefault (defParams->gradient.centerX);
	centerY->setDefault (defParams->gradient.centerY);

	if (pedited) {
		degree->setDefaultEditedState (pedited->gradient.degree ? Edited : UnEdited);
		feather->setDefaultEditedState (pedited->gradient.feather ? Edited : UnEdited);
		strength->setDefaultEditedState (pedited->gradient.strength ? Edited : UnEdited);
		centerX->setDefaultEditedState (pedited->gradient.centerX ? Edited : UnEdited);
		centerY->setDefaultEditedState (pedited->gradient.centerY ? Edited : UnEdited);
	} else {
		degree->setDefaultEditedState (Irrelevant);
		feather->setDefaultEditedState (Irrelevant);
		strength->setDefaultEditedState (Irrelevant);
		centerX->setDefaultEditedState (Irrelevant);
		centerY->setDefaultEditedState (Irrelevant);
	}
}

void Gradient::adjusterChanged (Adjuster* a, double newval) {

	updateGeometry (int(centerX->getValue()), int(centerY->getValue()), feather->getValue(), degree->getValue());

	if (listener && enabled->get_active()) {

		if (a == degree)
			listener->panelChanged (EvGradientDegree, degree->getTextValue());
		else if (a == feather)
			listener->panelChanged (EvGradientFeather, feather->getTextValue());
		else if (a == strength)
			listener->panelChanged (EvGradientStrength, strength->getTextValue());
		else if (a == centerX || a == centerY)
			listener->panelChanged (EvGradientCenter, Glib::ustring::compose ("X=%1\nY=%2", centerX->getTextValue(), centerY->getTextValue()));
	}
}

void Gradient::enabledChanged () {

	if (batchMode) {
		if (enabled->get_inconsistent()) {
			enabled->set_inconsistent (false);
			enaConn.block (true);
			enabled->set_active (false);
			enaConn.block (false);
		}
	else if (lastEnabled)
		enabled->set_inconsistent (true);

		lastEnabled = enabled->get_active ();
	}
	if (listener) {
		if (enabled->get_active())
			listener->panelChanged (EvGradientEnabled, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvGradientEnabled, M("GENERAL_DISABLED"));
	}
}

void Gradient::setAdjusterBehavior (bool degreeadd, bool featheradd, bool strengthadd, bool centeradd)
{
	degree->setAddMode(degreeadd);
	feather->setAddMode(featheradd);
	strength->setAddMode(strengthadd);
	centerX->setAddMode(centeradd);
	centerY->setAddMode(centeradd);
}

void Gradient::trimValues (rtengine::procparams::ProcParams* pp)
{
	degree->trimValue(pp->gradient.degree);
	feather->trimValue(pp->gradient.feather);
	strength->trimValue(pp->gradient.strength);
	centerX->trimValue(pp->gradient.centerX);
	centerY->trimValue(pp->gradient.centerY);
}

void Gradient::setBatchMode (bool batchMode)
{
	ToolPanel::setBatchMode (batchMode);
	degree->showEditedCB ();
	feather->showEditedCB ();
	strength->showEditedCB ();
	centerX->showEditedCB ();
	centerY->showEditedCB ();
}

void Gradient::setEditProvider (EditDataProvider* provider) {
	EditSubscriber::setEditProvider(provider);
}

void Gradient::editToggled () {
	if (edit->get_active()) {
		subscribe();
	}
	else
		unsubscribe();
}

CursorShape Gradient::getCursor(int objectID) {
	switch (objectID) {
	case (0):
	case (1):
		return CSMoveRotate;
	case (2):
	case (3):
		{
		int angle = degree->getIntValue();
		if (angle<-135 || (angle>=-45 && angle<=45) || angle>135)
			return CSMove1DV;
		return CSMove1DH;
		}
	case (4):
		return CSMove2D;
	default:
		return CSOpenHand;
	}
}

bool Gradient::mouseOver(int modifierKey) {
	EditDataProvider* editProvider = getEditProvider();
	if (editProvider && editProvider->object!=lastObject) {
		if (lastObject > -1) {
			if (lastObject == 2 || lastObject == 3) {
				EditSubscriber::visibleGeometry.at(2)->state = Geometry::NORMAL;
				EditSubscriber::visibleGeometry.at(3)->state = Geometry::NORMAL;
			}
			else
				EditSubscriber::visibleGeometry.at(lastObject)->state = Geometry::NORMAL;
		}
		if (editProvider->object > -1) {
			if (editProvider->object == 2 || editProvider->object == 3) {
				EditSubscriber::visibleGeometry.at(2)->state = Geometry::PRELIGHT;
				EditSubscriber::visibleGeometry.at(3)->state = Geometry::PRELIGHT;
			}
			else
				EditSubscriber::visibleGeometry.at(editProvider->object)->state = Geometry::PRELIGHT;
		}
		lastObject = editProvider->object;
		return true;
	}
	return false;
}

bool Gradient::button1Pressed(int modifierKey) {
	if (!(modifierKey & (GDK_CONTROL_MASK | GDK_CONTROL_MASK))) {
		// button press is valid (no modifier key)
		PolarCoord pCoord;
		EditDataProvider *provider = getEditProvider();
		int imW, imH;
		provider->getImageSize(imW, imH);
		double halfSizeW = imW/2.;
		double halfSizeH = imH/2.;
		draggedCenter.set(int(halfSizeW+halfSizeW*(centerX->getValue()/100.)), int(halfSizeH+halfSizeH*(centerY->getValue()/100.)));

		// trick to get the correct angle (clockwise/counter-clockwise)
		Coord p1 = draggedCenter;
		Coord p2 = provider->posImage;
		int p = p1.y;
		p1.y = p2.y;
		p2.y = p;

		pCoord.setFromCartesian(p1, p2);
		draggedPointOldAngle = pCoord.angle;
		//printf("\ndraggedPointOldAngle=%.3f\n\n", draggedPointOldAngle);
		draggedPointAdjusterAngle = degree->getValue();
		if (lastObject==2 || lastObject==3) {
			// Dragging a line to change the angle
			PolarCoord draggedPoint;
			Coord currPos;
			currPos = provider->posImage;
			Coord centerPos = draggedCenter;

			double diagonal = sqrt(double(imW)*double(imW)+double(imH)*double(imH));

			// trick to get the correct angle (clockwise/counter-clockwise)
			int p = centerPos.y;
			centerPos.y = currPos.y;
			currPos.y = p;

			draggedPoint.setFromCartesian(centerPos, currPos);
			// compute the projected value of the dragged point
			draggedFeatherOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue())/180.*M_PI);
			if (lastObject==3)
				draggedFeatherOffset = -draggedFeatherOffset;
			draggedFeatherOffset -= (feather->getValue() / 200. * diagonal);
		}
		return false;
	}
	else {
		// this will let this class ignore further drag events
		if (lastObject > -1) { // should theoretically always be true
			if (lastObject == 2 || lastObject == 3) {
				EditSubscriber::visibleGeometry.at(2)->state = Geometry::NORMAL;
				EditSubscriber::visibleGeometry.at(3)->state = Geometry::NORMAL;
			}
			else
				EditSubscriber::visibleGeometry.at(lastObject)->state = Geometry::NORMAL;
		}
		lastObject = -1;
		return true;
	}
}

bool Gradient::button1Released() {
	draggedPointOldAngle = -1000.;
	return true;
}

bool Gradient::drag(int modifierKey) {
	// compute the polar coordinate of the mouse position
	EditDataProvider *provider = getEditProvider();
	int imW, imH;
	provider->getImageSize(imW, imH);
	double halfSizeW = imW/2.;
	double halfSizeH = imH/2.;

	if (lastObject==0 || lastObject==1) {

		// Dragging a line to change the angle
		PolarCoord draggedPoint;
		Coord currPos;
		currPos = provider->posImage+provider->deltaImage;
		Coord centerPos = draggedCenter;

		// trick to get the correct angle (clockwise/counter-clockwise)
		int p = centerPos.y;
		centerPos.y = currPos.y;
		currPos.y = p;

		draggedPoint.setFromCartesian(centerPos, currPos);
		double deltaAngle = draggedPoint.angle - draggedPointOldAngle;
		if (deltaAngle>180.) // crossing the boundary (0->360)
			deltaAngle -= 360.;
		else if (deltaAngle<-180.) // crossing the boundary (360->0)
			deltaAngle += 360.;
		draggedPointOldAngle = draggedPoint.angle;

		draggedPointAdjusterAngle += deltaAngle;
		if (draggedPointAdjusterAngle > 180.)
			draggedPointAdjusterAngle = -360. + draggedPointAdjusterAngle;
		else if (draggedPointAdjusterAngle < -180.)
			draggedPointAdjusterAngle = 360. - draggedPointAdjusterAngle;
		//printf("draggedPointOldAngle: %.3f /  From %d,%d to %d,%d -> angle = %.3f  /  ", draggedPointAdjusterAngle, centerPos.x, centerPos.y, currPos.x, currPos.y, draggedPoint.angle);
		//printf("currAngle: %.3f = degree: %.3f + deltaAngle: %.3f %s /  draggedPointOldAngle: %.3f\n", draggedPointAdjusterAngle, degree->getValue(), deltaAngle, degree->getValue()>180.?">180":degree->getValue()<180.?"<180":"", draggedPointOldAngle);
		if (int(draggedPointAdjusterAngle) != degree->getIntValue()) {
			degree->setValue(draggedPointAdjusterAngle);
			updateGeometry (int(centerX->getValue()), int(centerY->getValue()), feather->getValue(), degree->getValue());
			if (listener)
				listener->panelChanged (EvGradientDegree, degree->getTextValue());
			return true;
		}
	}
	else if (lastObject==2 || lastObject==3) {
		// Dragging the upper or lower feather bar
		PolarCoord draggedPoint;
		Coord currPos;
		currPos = provider->posImage+provider->deltaImage;
		Coord centerPos = draggedCenter;

		double diagonal = sqrt(double(imW)*double(imW)+double(imH)*double(imH));

		// trick to get the correct angle (clockwise/counter-clockwise)
		int p = centerPos.y;
		centerPos.y = currPos.y;
		currPos.y = p;

		draggedPoint.setFromCartesian(centerPos, currPos);
		double currDraggedFeatherOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue())/180.*M_PI);
		if (lastObject==2)
			// Dragging the upper feather bar
			currDraggedFeatherOffset -= draggedFeatherOffset;
		else if (lastObject==3)
			// Dragging the lower feather bar
			currDraggedFeatherOffset = -currDraggedFeatherOffset + draggedFeatherOffset;
		currDraggedFeatherOffset = currDraggedFeatherOffset * 200. / diagonal;

		if (int(currDraggedFeatherOffset) != feather->getIntValue()) {
			feather->setValue(double(int(currDraggedFeatherOffset)));
			updateGeometry (centerX->getValue(), centerY->getValue(), feather->getValue(), degree->getValue());
			if (listener)
				listener->panelChanged (EvGradientFeather, feather->getTextValue());
			return true;
		}
	}
	else if (lastObject==4) {
		// Dragging the circle to change the center
		Coord currPos;
		draggedCenter += provider->deltaPrevImage;
		currPos = draggedCenter;
		currPos.clip(imW, imH);
		int newCenterX = int((double(currPos.x)-halfSizeW)/halfSizeW*100.);
		int newCenterY = int((double(currPos.y)-halfSizeH)/halfSizeH*100.);
		if (newCenterX!=centerX->getIntValue() || newCenterY!=centerY->getIntValue()) {
			centerX->setValue(newCenterX);
			centerY->setValue(newCenterY);
			updateGeometry (newCenterX, newCenterY, feather->getValue(), degree->getValue());
			if (listener)
				listener->panelChanged (EvGradientCenter, Glib::ustring::compose ("X=%1\nY=%2", centerX->getTextValue(), centerY->getTextValue()));
			return true;
		}
	}
	return false;
}

void Gradient::switchOffEditMode () {
	if (edit->get_active()) {
		// switching off the toggle button
		bool wasBlocked = editConn.block(true);
		edit->set_active(false);
		if (!wasBlocked) editConn.block(false);
	}
	EditSubscriber::switchOffEditMode();  // disconnect
}

