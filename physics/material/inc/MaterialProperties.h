
#ifndef MATERIALPROPERTIES_H
#define MATERIALPROPERTIES_H

#include <iostream>

namespace geant{
  /**
   * @brief   Class to represent derived material properties.
   * @class   MaterialProperties
   * @author  M Novak, A Ribon
   * @date    January 2016
   *
   * Each material has a corresponding material properties object that initialised
   * based on the material. The material properties object stores information like
   * total number of atoms per volume, total number of electrons per volume,
   * number of atoms per volume for each constituting element, mean excitation
   * energy and parameters that are needed to compute the density effect parameter.
   *
   * The mean excitation energy (\f$ I \f$ ):
   *  - is taken from the internal NISTMaterialData database if the material
   *    can be identified in the internal NISTMaterialData database based on its
   *    name and the mean excitation energy stored in the database is non-zero
   *  - otherwise, it is computed from the elemental mean excitation energies
   *    (taken from the NISTMaterialData database) by using the relation
   *    \f[
   *       \ln(I) = \frac{\sum_i p_i Z_i \ln(I(Z_i))}{\sum_i p_i Z_i}
   *    \f]
   *    where \f$ p_i \f$ is the proportion by number of the \f$ i\f$-th element
   *    in the material (i.e. number of \f$ i\f$-th atoms per volume) with atomic
   *    number \f$ Z_i \f$ and mean excitation energy \f$ I(Z_i) \f$.
   *
   * The density effect parameter \f$ \delta(X) \f$ where \f$ X = \log_{10}(pc/(m_0 c^2)) = \log_{10}(\gamma\beta) \f$
   * (with \f$ pc \f$ total momentum in energy units, \f$ m_0 c^2 \f$ rest mass energy,
   *  \f$ \beta=v/c \f$ and \f$ \gamma = (1-\beta^2)^{-1/2} \f$):
   * -# is computed according to the parametrisation given in \cite sternheimer1984density (Eqs.(14,9,10))
   *    if the material can be identified in the internal DensityEffectData database based on its
   *    name (see \ref DensityEffectDataDoc). Namely
   *    \f[
   *     \delta(X) =
   *       \begin{cases}
   *        \delta_0 = \delta(X_0) 10^{2(X-X_0)}   & \quad  X \leq X_0   \\
   *        2\ln(10)X + a(X_1-X)^m - \tilde{C}     & \quad  X_0 < X <X_1 \\
   *        2\ln(10)X - \tilde{C}                  & \quad  X \geq X_1
   *       \end{cases}
   *    \f]
   *    where \f$ \delta(X_0), X_0, X_1, a, m \f$ parameters are taken from the DensityEffectData database,
   *    \f$\tilde{C} = -C = 2\ln(I/(h\nu_p)) + 1\f$ where \f$ I \f$ is the mean excitation energy (see above)
   *    and \f$ h\nu_{p} \f$ is the plasma energy of the electrons of the material considered as free
   *    electron gas. Parameters \f$ \tilde{C} = -C, h\nu_{p} \f$ are also taken from the DensityEffectData
   *    database if the material can be identified based on its name. Note, that the density effect parameters
   *    were determined in \cite sternheimer1984density through fitting the above formula to computed \f$ \delta(X)\f$
   *    values for the given materials.
   * -# if the material cannot be identified in the DensityEffectData database based on its name:
   *    1. plasma energy is computed
   *      \f[
   *         h\nu_{p} = \sqrt{4\pi n_{el} r_e} \hbar c
   *      \f]
   *      where \f$ n_{el}\f$ is the electron density of the material and \f$ r_e \f$ is the classical electron radius.
   *      The electron density for simple (elemental) materials
   *      \f[
   *         n_{el} = Z \frac{\mathcal{N}_{av}\rho}{A}
   *      \f]
   *      and compound materials
   *      \f[
   *         n_{el} = \sum_{i}^{\#elements} Z_i \frac{\mathcal{N}_{av}w_i\rho}{A_i}
   *      \f]
   *      where \f$ \mathcal{N}_{av} \f$ is the Avogadro number, \f$ \rho \f$ is the material density,
   *      \f$ w_i \f$ is the fractional mass of the \f$ i \f$-th element of the material
   *      with molar mass of \f$ A_i \f$.
   *    2. then the general formula given in \cite sternheimer1971general is used to compute the density
   *       effect parameters. Namely
   *      \f[
   *       \delta(X) =
   *         \begin{cases}
   *          \delta(X_0) = 0                        & \quad  X \leq X_0   \\
   *          2\ln(10)X + a(X_1-X)^m - \tilde{C}     & \quad  X_0 < X <X_1 \\
   *          2\ln(10)X - \tilde{C}                  & \quad  X \geq X_1
   *         \end{cases}
   *      \f]
   *      where \f$\tilde{C} = -C = 2\ln(I/(h\nu_p)) + 1\f$ as mentioned above, \f$ m = 3\f$ is fixed.
   *      The value of \f$ a \f$ is determined by the requirement
   *      \f$ \delta(X=X_0) = 0 = 2\ln(10)X_0 + a(X_1-X_0)^m - \tilde{C} \f$ and gives
   *      \f$ a = [\tilde{C}-2\ln(10)X_0]/[X_1-X_0]^m\f$ that can be written as
   *      \f$ a = 2\ln(10)[X_a-X_0]/[X_1-X_0]^m\f$ with \f$ X_a=\tilde{C}/(2\ln(10))\f$.
   *      The values of the remaining parameters \f$ X_0, X_1 \f$ according to \cite sternheimer1971general
   *      - for solids and liquids(Eqs.(33-36) in \cite sternheimer1971general)
   *      \f[
   *         I < 100 \textrm{[eV]}
   *         \begin{cases}
   *          X_0 = 0.2,                 X_1 = 2.0   & \textrm{for} \quad  \tilde{C}  <   3.681  \\
   *          X_0 = 0.326\tilde{C}-1.0,  X_1 = 2.0   & \textrm{for} \quad  \tilde{C} \geq 3.681  \\
   *         \end{cases}
   *      \f]
   *      \f[
   *         I \geq 100 \textrm{[eV]}
   *         \begin{cases}
   *          X_0 = 0.2,                 X_1 = 3.0   & \textrm{for} \quad  \tilde{C}  <   5.215  \\
   *          X_0 = 0.326\tilde{C}-1.5,  X_1 = 3.0   & \textrm{for} \quad  \tilde{C} \geq 5.215  \\
   *         \end{cases}
   *      \f]
   *      - for gases(Eqs.(37-41) plus text in \cite sternheimer1971general)
   *      \f[
   *         \begin{matrix}
   *          X_0 &=& 1.6,                 & X_1 = 4.0   & \textrm{for} &        &      & \tilde{C} & < & 10.0  \\
   *          X_0 &=& 1.7,                 & X_1 = 4.0   & \textrm{for} & 10.0   & \leq & \tilde{C} & < & 10.5  \\
   *          X_0 &=& 1.8,                 & X_1 = 4.0   & \textrm{for} & 10.5   & \leq & \tilde{C} & < & 11.0  \\
   *          X_0 &=& 1.9,                 & X_1 = 4.0   & \textrm{for} & 11.0   & \leq & \tilde{C} & < & 11.5  \\
   *          X_0 &=& 2.0,                 & X_1 = 4.0   & \textrm{for} & 11.5   & \leq & \tilde{C} & < & 12.25  \\
   *          X_0 &=& 2.0,                 & X_1 = 5.0   & \textrm{for} & 12.25  & \leq & \tilde{C} & < & 13.804  \\
   *          X_0 &=& 0.326\tilde{C}-2.5,  & X_1 = 5.0   & \textrm{for} & 13.804 & \leq & \tilde{C} &   &
   *         \end{matrix}
   *      \f]
   * -# the last point is to account the fact that the density effect parameters for gaseous materials were derived under
   *    the conditions \f$ T = 0[^{\circ}\textrm{C}]\;\textrm{and}\; P = 1 \textrm{[atm]}\f$ and the corresponding gas
   *    density is \f$ \rho_0 \f$. According to Eq.(44) in \cite sternheimer1971general, if the gas density is
   *    \f$ \rho=\eta\rho_0 \f$, we need to evaluate \f$ \delta \f$ at \f$ p\eta^{1/2}\f$ momentum instead of \f$ p \f$.
   *    Since \f$ X=\log_{10}(p/(m_0 c)) = \log_{10}(p)-\log_{10}(m_0c)\f$
   *    \f$ \to X_{\eta} \equiv \log_{10}(p)+0.5\log_{10}(\eta)-\log_{10}(m_0c) = X+0.5\log_{10}(\eta)\f$, elavulating
   *    \f$ \delta \f$ at \f$ p\eta^{1/2}\f$ momentum corresponds to \f$ \delta(X_{\eta}\equiv X+0.5\log_{10}(\eta)) \f$.
   *    Therefore, evaluating \f$ \delta \f$ at momentum \f$ p\eta^{1/2}\f$ instead of \f$ p \f$ corresponds to taking
   *    the value of \f$ \delta \f$ at a shifted \f$ X \f$. Since \f$ X \f$ is shifted to \f$ X+0.5\log_{10}(\eta) \f$
   *    it corresponds to the shift of the two limits in \f$ X \f$ i.e.
   *      \f[
   *         \begin{matrix}
   *              X_{0} & \to & X_{0,\eta}       & \equiv & X_0 - 0.5\log_{10}(\eta) \\
   *              X_{1} & \to & X_{1,\eta}       & \equiv & X_1 - 0.5\log_{10}(\eta)
   *         \end{matrix}
   *      \f]
   *    and
   *    \f[
   *       \delta(X_{\eta}\equiv X+0.5\log_{10}(\eta)| X_0,X_1,\tilde{C},a,m)
   *       = 2\ln(10)X+\ln(10)\log_{10}(\eta) - \tilde{C} + a [X_1-X-0.5\log_{10}(\eta)]^m
   *       = 2\ln(10)X- \tilde{C}_{\eta} + a [X_{1,\eta}-X]^m
   *       = \delta(X| X_{0,\eta},X_{1,\eta},\tilde{C}_{\eta},a,m)
   *    \f]
   *    where
   *      \f[
   *              \tilde{C}_{\eta} \equiv  \tilde{C} - \ln(10)\log_{10}(\eta) \\
   *      \f]
   *    It can be seen from the condition for \f$ a \f$ that parameters \f$ a,\; m \f$ remains the same i.e
   *      \f[
   *        a = \frac{\tilde{C}-2\ln(10)X_0}{[X_1-X_0]^m}
   *          = \frac{\tilde{C}_{\eta}-2\ln(10)X_{0,\eta}}{[X_{1,\eta}-X_{0,\eta}]^m}
   *      \f]
   *    and \f$ X_a \f$ now becomes \f$ X_{a,\eta} \equiv \tilde{C}_{\eta}/(2\ln(10)) \f$.\\
   *    It is important to note, that this correction can also be applied to solids and liquids that has a
   *    density different from that was used to derive the density effect parameters (ro data).
   *
   *\todo
   * -# Add possibility to set mean excitation energy.
   * -# Add possibility to set density effect parameters directly (X0,X1,a,m,Cbar).
   * -# The best would be to use numerical \f$ \delta(X) \f$ values that were used to derive the fit parameters.
   *    Such data avaliable for a significant set of simple and composite materials (EGS simulations use this)
   *    If we take those data, the data file will contain the mean ionisation energy (I[eV]), the density (\f$ \rho_0 \f$)
   *    that were used in the computation and the \f$ E_{kin},\; \delta(E_{kin}) \f$ datapairs for electrons. What we need
   *    to do is:
   *    - first we need to change the \f$ E_{kin} \f$ scale to \f$ X=log_{10}(pc/(m_0c^2))\f$. Since
   *      \f$ (pc/(m_0c^2)) = \sqrt{ E_{kin}/(m_0c^2)[E_{kin}/(m_0c^2)+2]} \f$ where \f$ m_0 c^2 \f$ is the electron rest mass energy,
   *      \f$ X = log_{10} (\sqrt{ E_{kin}/(m_0c^2)[E_{kin}/(m_0c^2)+2]}) \f$.
   *    - then we need to apply possible correction comming from the possibility that the density hat was assumed
   *      during the computation of these \f$ \delta(E_{kin}) \f$ data (\f$ \rho_0 \f$) and our current density (\f$ \rho \f$)
   *      are different. We can do this by simple shifting the \f$ X \f$ scale computed above to \f$ X = X - 0.5\log_{10}(\eta) \f$
   *      where \f$ \eta = \rho/\rho_0 \f$.
   *    - it's also interesting to note, that in the upper limiting part of \f$ X \f$, \f$ \delta(X) = 2\ln(10)X - \tilde{C}\f$
   *      where \f$\tilde{C} = 2\ln(I/(h\nu_p)) + 1\f$. Plugging this later into the previous gives
   *      \f$ I/(h\nu_p) = \exp{[X\ln(10)-0.5(\delta(X)+1)]} \f$ i.e. taking the \f$ X,\; \delta(X) \f$  pair that corresponds
   *      to the highest value of \f$ X \f$ gives the possiblity of getting \f$ h\nu_p \f$ if the \f$ I \f$ was given in the file.
   *
   */
class Material;

class MaterialProperties{
 public:
/**
* @name Constructors/destructors:
*/
//@{
   /** @brief Constructor to create material properties object.
    *
    *  Initialisation of parameters are done at construction time.
    *
    *  @param[in] Pointer to the corresponding material object.
    */
   MaterialProperties(Material *mat);
   /** @brief Destructor.
    */
  ~MaterialProperties();
//@}

/**
* @name Basic material property getters:
*/
//@{
    /** @brief Public method to get the total number of atoms per volume in this material.
     *  @return Total  number of atoms per volume in this material.
     */
    double         GetTotalNumOfAtomsPerVol() const {return fTotalNumOfAtomsPerVol;}

    /** @brief Public method to get the total number of aelectrons per volume in this material.
     *  @return Total number of electrons per volume in this material.
     */
    double         GetTotalNumOfElectronsPerVol() const {return fTotalNumOfElectronsPerVol;}

    /** @brief Public method to get the number of atoms per volume for each elements this material is built up.
     *  @return Number of atoms per volume for each elements this material is built up [Material::GetNumberOfElements()]. */
    const double*   GetNumOfAtomsPerVolumeVect() const {return fNumOfAtomsPerVolVect;}
//@}


/**
* @name Density effect parameter:
*/
//@{
      /**
       * @brief Public method to compute density effect correction parameter \f$ \delta(X)\f$.
       * @param[in]  val This is \f$ X\equiv \log_{10} (\gamma\beta) = \ln(\gamma^2\beta^2)/(2\ln(10))\f$.
       * @return     Density effect correction parameter \f$ \delta(X)\f$.
       */
      double GetDensityEffectFunctionValue(const double val);
//@}


/**
 * @name Ionization parameter getters/setters:
 */
//@{
   /** @brief Public method to get the mean excitation energy of this material.
    *  @return Mean excitation energy in internal [energy] units.              */
   double  GetMeanExcitationEnergy() const { return fMeanExcitationEnergy;}
//@}

/**
 * @name Others:
 */
//@{
   /** @brief Public method to get the plasma energy in internal units.
    *  @return Plasma energy in internal [energy] units.              */
   double  GetPlasmaEnergy() const { return fPlasmaEnergy;}

   /** @brief Public method to get the radiation length in internal units.
    *
    *  Radiation length is computed at initialisation based on \cite tsai1974pair
    *  (Eqs.(3.65) with Coulomb correction from \cite davies1954theory Eqs.(36-38).
    *
    *  @return Radiation length in internal [length] units.              */
   double  GetRadiationLength() const { return fRadiationLength;}
//@}

/**
 * @name Printouts:
 */
//@{
   friend std::ostream& operator<<(std::ostream&, const MaterialProperties*);
   friend std::ostream& operator<<(std::ostream&, const MaterialProperties&);
//@}


 private:
   void InitialiseMembers();
   /** @brief Internal method to compute some frequently used basic material properties.*/
   void ComputeBasicMaterialParameters();
   void ComputeIonizationParameters();
   void ComputeDensityEffectParameters();
   void ComputeRadiationLength();

 private:
   /** @brief Pointer to the material that this material parameters belong to */
   Material *fMaterial;

/**
* @name Basic material parameters:
*/
//@{
   /** @brief Number of atoms per volume for each element this material is buit up. */
   double       *fNumOfAtomsPerVolVect;
   /** @brief Total number of atoms per volume in this material. */
   double        fTotalNumOfAtomsPerVol;
   /** @brief Total number of electrons per volume in this material. */
   double        fTotalNumOfElectronsPerVol;
//@}
/**
* @name Ionization parameter:
*/
//@{
   /** @brief Mean excitation energy of this material in internal [energy] units. */
   double        fMeanExcitationEnergy;
//@}
//@}
/**
* @name Density effect parameters:
*/
//@{
   /** @brief Plasma energy in internal units. */
   double fPlasmaEnergy;
   /** @brief Density effect parameter Cbar. */
   double fDensityEffectParameterC;
   /** @brief Density effect parameter X_0. */
   double fDensityEffectParameterX0;
   /** @brief Density effect parameter X_1. */
   double fDensityEffectParameterX1;
   /** @brief Density effect parameter a. */
   double fDensityEffectParameterA;
   /** @brief Density effect parameter m. */
   double fDensityEffectParameterM;
   /** @brief Density effect parameter delta_0. */
   double fDensityEffectParameterDelta0;
   /** @brief Radiation length in internal [length] units */
   double fRadiationLength;
//@}


};

}    // namespace geant

#endif // MATERIALPROPERTIES_H
