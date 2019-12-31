#include "Jet.h"

bool passPUID(analysis::core::Jet j, int year)
{
   bool loosePU = bool(j._fullid & (1 << 2)) or j._pt > 50;
   bool mediumPU = bool(j._fullid & (1 << 1)) or j._pt > 50;
   bool tightPU = bool(j._fullid & (1 << 0));

   if (year == 2016 or year == 2018)
      return loosePU;
   else if (year == 2017)
   {
      if (TMath::Abs(j._eta) >= 2.60 and TMath::Abs(j._eta) <= 3.0)
      {
         return tightPU;
      }
      return loosePU;
   }
   // shouldn't happen
   std::cout << "year invalid" << std::endl;    
   return false; 
}


bool passJetID(analysis::core::Jet j, int year)
{
   // for 2016 use Loose, for 2017 and 2018 use tight.. rgerosa 15/11
   bool tightID = false;
   bool looseID = false;

   double jeta = TMath::Abs(j._eta);
   int numConst = j._cm + j._nm;

   if (year == 2016)
   {
      if (jeta <= 2.7)
      {
         looseID = (j._nhf < 0.99 and j._nef < 0.99 and numConst > 1);
         // tightID = (j._nhf < 0.90 and j._nef < 0.90 and numConst > 1);

         if (jeta <= 2.4)
         {
            // tightID &= (j._chf > 0 and j._cm > 0 and j._cef < .99);
            looseID &= (j._chf > 0 and j._cm > 0 and j._cef < .99);
         }
      }
      else if (jeta <= 3.0)
      {
         // tightID = (j._nef > 0.01 and j._nhf < 0.98 and j._nm > 2);
         looseID = (j._nef > 0.01 and j._nhf < 0.98 and j._nm > 2);
      }
      else
      {
         // tightID = (j._nef < 0.90 and j._nm > 10);
         looseID = (j._nef < 0.90 and j._nm > 10);
      }

      return looseID;
   }

   else if (year == 2017)
   {
      if (jeta <= 2.7)
      {
         tightID = (j._nhf < 0.90 and j._nef < 0.90 and numConst > 1);

         if (jeta <= 2.4)
         {
            tightID &= (j._chf > 0 and j._cm > 0);
         }
      }
      else if (jeta <= 3.0)
      {
         tightID = (j._nef > 0.02 and j._nef < 0.99 and j._nm > 2);
      }
      else
      {
         tightID = (j._nef < 0.90 and j._nhf > 0.02 and j._nm > 10);
      }

      return tightID;
   }

   else if (year == 2018)
   {
      if (jeta <= 2.6)
      {
         tightID = (j._nhf < 0.90 and j._nef < 0.90 and numConst > 1 and j._chf > 0 and j._cm > 0 );
         // and j._muf < 0.80 and j._cef < 0.80 lepVeto recommended
      }

      else if (jeta <= 2.7) 
      {
         tightID = (j._nhf < 0.90 and j._nef < 0.99  and j._cm > 0 );
         // and j._cef < 0.80 and j._muf < 0.80 lepVeto recommended
      }

      else if (jeta <= 3.0)
      {
         tightID = (j._nef > 0.02 and j._nef < 0.99 and j._nm > 2);
      }
      else
      {
         tightID = (j._nhf > 0.2 and j._nef < 0.90 and j._nm > 10);

      }
      return tightID;
   }   
   std::cout << "year invalid" << std::endl;    
   return false; 
}

// no longer needed?
// bool passNoiseJet(analysis::core::Jet j)
// {
//    if (!mcLabel and year == 2017)
//    {
//       float jeta = TMath::Abs(j._eta);
//       float jpt = j._pt;
//       if (jeta >= 2.65 and jeta <= 3.139 and jpt < 50)
//          return false;
//    }
//    return true;
// }

// try to use the bit from the fullId, this is valid for 8X training only
// bool passLoosePUID(analysis::core::Jet j)
// {
//    float jeta = TMath::Abs(j._eta);
//    float jpt = j._pt;
//    float jpuid = j._puid;

//    if (jeta < 2.5)
//    {
//       if (jpt >= 30 and jpt < 50 and jpuid < -0.89)
//          return false;
//       if (jpt >= 20 and jpt < 30 and jpuid < -0.97)
//          return false;
//    }
//    else if (jeta < 2.75)
//    {
//       if (jpt >= 30 and jpt < 50 and jpuid < -0.52)
//          return false;
//       if (jpt >= 20 and jpt < 30 and jpuid < -0.68)
//          return false;
//    }
//    else if (jeta < 3.0)
//    {
//       if (jpt >= 30 and jpt < 50 and jpuid < -0.38)
//          return false;
//       if (jpt >= 20 and jpt < 30 and jpuid < -0.53)
//          return false;
//    }
//    else if (jeta < 5)
//    {
//       if (jpt >= 30 and jpt < 50 and jpuid < -0.30)
//          return false;
//       if (jpt >= 20 and jpt < 30 and jpuid < -0.47)
//          return false;
//    }
//    return true;
// }

