//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "HadesTestApp.h"
#include "HadesApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
HadesTestApp::validParams()
{
  InputParameters params = HadesApp::validParams();
  return params;
}

HadesTestApp::HadesTestApp(InputParameters parameters) : MooseApp(parameters)
{
  HadesTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

HadesTestApp::~HadesTestApp() {}

void
HadesTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  HadesApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"HadesTestApp"});
    Registry::registerActionsTo(af, {"HadesTestApp"});
  }
}

void
HadesTestApp::registerApps()
{
  registerApp(HadesApp);
  registerApp(HadesTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
HadesTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  HadesTestApp::registerAll(f, af, s);
}
extern "C" void
HadesTestApp__registerApps()
{
  HadesTestApp::registerApps();
}
