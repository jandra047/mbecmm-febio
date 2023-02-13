// This file defines the required functions for registering the plugin classes
// with the FEBio framework.

#include "stdafx.h"
#include "FECore/FECoreKernel.h"
#include "FEMbeCmm.h"

#ifdef WIN32
#include "targetver.h"
#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
	switch (ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
	return TRUE;
}
#endif // WIN32


//-----------------------------------------------------------------------------
// This template class instantiates a factory class that will be used by FEBio
// to create an instance of the plugin class. It is this factory class that needs
// to be registered with the FEBio framework. This is done by returing a pointer
// to this factory object in the RegisterPlugin function below.
FEPluginFactory_T<FEMbeCmm, FEMATERIAL_ID> mbe_cmm_factory("mbe_cmm");

//-----------------------------------------------------------------------------
// This required function returns the version of the FEBio SDK that is being
// used by the plugin. This version number will be checked by FEBio to make
// sure that the FEBio executable is compatible with this SDK. Usually this
// function just returns the predefined macro FE_SDK_VERSION.
FECORE_EXPORT unsigned int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT void GetPluginVersion(int& major, int& minor, int& patch)
{
	major = 1;
	minor = 33;
	patch = 7;
}

//-----------------------------------------------------------------------------
// If FEBio is successful in loading the plugin, this is the first function
// that will be called and can be used to initialize any resources that are
// needed by the plugin. FEBio passes a reference to the FECoreKernel as a parameter,
// which can be used to access some of the FEBio resources such as the log file.
// This function is optional.
FECORE_EXPORT void PluginInitialize(FECoreKernel& febio)
{
	FECoreKernel::SetInstance(&febio);
}

//-----------------------------------------------------------------------------
// This required function should return the number of classes that this plugin defines.
// This number determines the number of times that FEBio calls the RegisterPlugin
// later, which is used to register the actual factory classes.
FECORE_EXPORT int PluginNumClasses()
{
	return 1;
}

//-----------------------------------------------------------------------------
// This required function is called by FEBio and will be called the same number of times
// as return by PluginNumClasses. The parameter i is used as a counter which contains
// how often this function has been called thus far.
FECORE_EXPORT FECoreFactory* PluginGetFactory(int i)
{
	if (i == 0) return &mbe_cmm_factory;
	return 0;
}

//-----------------------------------------------------------------------------
// This function is called when FEBio exits and gives the plugin a chance to close
// any resource that have been allocated. 
// This is an optional function.
FECORE_EXPORT void PluginCleanup()
{

}