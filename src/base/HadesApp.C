#include "HadesApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
HadesApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  return params;
}

HadesApp::HadesApp(InputParameters parameters) : MooseApp(parameters)
{
  HadesApp::registerAll(_factory, _action_factory, _syntax);
}

HadesApp::~HadesApp() {}

void
HadesApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"HadesApp"});
  Registry::registerActionsTo(af, {"HadesApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
HadesApp::registerApps()
{
  registerApp(HadesApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
HadesApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  HadesApp::registerAll(f, af, s);
}
extern "C" void
HadesApp__registerApps()
{
  HadesApp::registerApps();
}
