#pragma once

#include "Interfaces.h"
#include "Models.h"
#include <cmath>
#include <map>
#include <memory>
#include <string>

namespace mc {

// Canonical key for the correlation between two risk factors, e.g.
// corrKey("EUR.USD.0", "GBP.USD.0") -> "CORR.EUR.USD.0:GBP.USD.0"
// The key is order independent (names are sorted), so setCorrelation(A,B) == setCorrelation(B,A).
// It deliberately does NOT contain the model family: one correlation serves every model of the asset.
inline std::string corrKey(const std::string& rfA, const std::string& rfB) {
    return rfA < rfB ? "CORR." + rfA + ":" + rfB : "CORR." + rfB + ":" + rfA;
}

// Model tag "<FAMILY>.<ASSET>" -> asset (e.g. "BS_EULER.EUR.USD" -> "EUR.USD")
// Everything after the first '.' is the asset, because asset names themselves contain a dot.
inline std::string assetOfModelTag(const std::string& tag) {
    auto dot = tag.find('.');
    MYASSERT(dot != std::string::npos, "Invalid model tag format: " << tag);
    return tag.substr(dot + 1);
}

// ----------------------------------------------------------------------------
// Market: everything the engine needs to know about the world on the valuation date.
//   - discount curves, one per currency (USD is the numeraire)
//   - FX spots & vols of the simulated "<CCY>.USD" assets
//   - correlations between risk factors ("<asset>.<brownian index>")
//   - model factory: builds the diffusion models from market data and a PricerConfig
// FX drifts are not inputs: they are implied by the two discount curves, which keeps
// the simulation consistent with the numeraire (E[S(T)] = forward).
// ----------------------------------------------------------------------------
// The Market object is a read-only snapshot used by the engine. It replaces the former
// "MarketDataServer": curves and models are now part of the market and not of a bag of maps.
// Invariants:
//   - all quantities refer to the valuation date `today`;
//   - a missing curve / spot / vol is an error (MYASSERT), never a silent default;
//   - only "<CCY>.USD" assets are quoted and simulated (USD is the numeraire); other pairs are
//     derived by the pricers (see PairObservable in Pricers.h).
class Market {
public:
    explicit Market(Date today = 0) : m_today(today) {}

    Date today() const { return m_today; }

    // --- Curves -------------------------------------------------------------
    void setCurve(const std::string& ccy, std::shared_ptr<const ICurve> curve) { m_curves[ccy] = std::move(curve); }
    // Convenience for demos and tests: installs a FlatDiscountCurve for the currency.
    void setFlatCurve(const std::string& ccy, double rate) {
        setCurve(ccy, std::make_shared<FlatDiscountCurve>("IR.DISCOUNT." + ccy, m_today, rate));
    }
    // Curve of a currency; fails loudly if absent.
    const ICurve& curve(const std::string& ccy) const {
        auto it = m_curves.find(ccy);
        MYASSERT(it != m_curves.end(), "Discount curve not found: " << ccy);
        return *it->second;
    }
    double df(const std::string& ccy, Date d) const { return curve(ccy).df(d); }

    // --- FX assets ("EUR.USD" = USD per EUR) --------------------------------
    void   setSpot(const std::string& asset, double s) { m_spots[asset] = s; }
    void   setVol(const std::string& asset, double v) { m_vols[asset] = v; }
    double spot(const std::string& asset) const { return lookup(m_spots, "Spot", asset); }
    double vol(const std::string& asset) const { return lookup(m_vols, "Vol", asset); }

    // Risk neutral drift of FX asset BASE.QUOTE in the QUOTE measure: r_quote - r_base,
    // read from the curves at a one year horizon (exact for flat curves).
    // Risk neutral drift of the FX rate BASE.QUOTE, expressed in the QUOTE measure:  r_quote - r_base.
    // Why: the rate is the QUOTE price of a BASE asset that earns r_base, so under the QUOTE
    // risk-neutral measure its expected growth is r_quote - r_base. This is implied by the two
    // curves, which keeps simulation, discounting and FX conversion mutually consistent (no arbitrage).
    // It is read from the curves as log(P_base(0,1y) / P_quote(0,1y)): exact for flat curves; a term
    // structure would need a time-dependent drift (extension).
    double fxDrift(const std::string& asset) const {
        auto dot = asset.find('.');
        MYASSERT(dot != std::string::npos, "Invalid FX asset: " << asset);
        const std::string base = asset.substr(0, dot), quote = asset.substr(dot + 1);
        const Date        h    = m_today + 365;
        return std::log(df(base, h) / df(quote, h));
    }

    // --- Correlations -------------------------------------------------------
    // Correlation between the first Brownian of two assets
    // Risk factors are named "<ASSET>.<index of the Brownian>" (e.g. "EUR.USD.0"). All models of
    // this example have ONE Brownian, so the correlation of two assets is the one between their
    // Brownians number 0; this helper hides the ".0".
    void setCorrelation(const std::string& assetA, const std::string& assetB, double rho) {
        m_correlations[corrKey(assetA + ".0", assetB + ".0")] = rho;
    }
    // Correlation between two risk factors used by the engine to build Sigma:
    // 1 on the diagonal, the stored value if present, 0 (independent) otherwise.
    double correlation(const std::string& rfA, const std::string& rfB) const {
        if (rfA == rfB)
            return 1.0;
        auto it = m_correlations.find(corrKey(rfA, rfB));
        return it != m_correlations.end() ? it->second : 0.0;
    }

    // --- Model factory ------------------------------------------------------
    // Model factory. `model_tag` = "<FAMILY>.<ASSET>"; market data (spot, vol, drift) comes from this
    // market, numerical parameters (e.g. the Euler step) from the config.
    // A new model is created on every call and not cached: models hold per-timeline state
    // (configureTimeline), and a cached one would be stale after a market bump (e.g. a Greek).
    std::shared_ptr<IModelMC> createModel(const std::string& model_tag, const PricerConfig& config) const {
        auto              dot    = model_tag.find('.');
        const std::string family = model_tag.substr(0, dot);
        const std::string asset  = assetOfModelTag(model_tag);
        const AssetConfig ac     = config.assetConfig(asset);

        const double S0 = spot(asset), mu = fxDrift(asset), sigma = vol(asset);
        switch (parseModelFamily(family)) {
        case ModelFamily::BS_CLOSED:
            return std::make_shared<ModelBSClosedForm>(model_tag, S0, mu, sigma);
        case ModelFamily::BS_EULER:
            return std::make_shared<ModelBSEuler>(model_tag, S0, mu, sigma, ac.euler_max_step);
        default:
            MYASSERT(false, "Unknown model family: " << family);
            return nullptr;
        }
    }

private:
    // Map lookup that reports WHAT is missing (spot, vol) and for which asset.
    static double lookup(const std::map<std::string, double>& m, const char* what, const std::string& key) {
        auto it = m.find(key);
        MYASSERT(it != m.end(), what << " not found in Market: " << key);
        return it->second;
    }

    Date                                           m_today;
    std::map<std::string, std::shared_ptr<const ICurve>> m_curves;
    std::map<std::string, double>                  m_spots, m_vols, m_correlations;
};

} // namespace mc
