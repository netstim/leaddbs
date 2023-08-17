'''
    By K. Butenko
    Functions for PathwayTune (see description in the headers)
'''


import json
import os
import sys

def create_NB_dictionaries(side, FF_dictionary, disease='spontaneous human combustion'):

    ''' Either imports Activation Profile Dictionary from Fiber Filtering (FF_dictionary)
        or takes a pre-defined dictionary for the given disease from TractSymptomLibrary
        and saves them in separate dictionaries for modulated symptoms, "soft" side-effects and "hard" side-effects
        in the stimulation folder / NB_ '''

    # create an output folder in the stim folder of the patient
    try:
        os.makedirs(os.environ['STIMDIR'] + '/NB_' + str(side))
    except:
        print("NB folder already exists")

    # retrieve activation profiles for the specific disease
    # "side" in the name does not actually play the role, it will be checked "on site"
    if FF_dictionary != 'no dictionary':
        with open(FF_dictionary, 'r') as fp:
            profiles = json.load(fp)
        fp.close()

        if 'profile_dict' in profiles.keys():
            with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/profile_dict.json', 'w') as save_as_dict:
                json.dump(profiles['profile_dict'], save_as_dict)
            profile_dict = profiles['profile_dict']
        else:
            print("profile_dict is missing, exiting...")
            raise SystemExit

        if 'Soft_SE_dict' in profiles.keys():
            Soft_SE_dict = profiles['Soft_SE_dict']
        else:
            print("Soft_SE_dict is missing")
            print("Will continue, but consider exporting negative tracts from FF")
            Soft_SE_dict = {}

        with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Soft_SE_dict.json', 'w') as save_as_dict:
            json.dump(Soft_SE_dict, save_as_dict)

        # IMPORTANT: What about SE_dict? Take it from some prior studies?
        if 'SE_dict' in profiles.keys():
            SE_dict = profiles['SE_dict']
        else:
            print("SE_dict is missing")
            print("Will continue, but consider exporting negative tracts from FF")
            SE_dict = {}

        with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/SE_dict.json', 'w') as save_as_dict:
            json.dump(SE_dict, save_as_dict)

    else:
        from TractSymptomLibrary import get_disease_profiles
        [profile_dict, Soft_SE_dict, SE_dict] = get_disease_profiles(disease)
        # here we need to implement a selection mechanism (but no adjustment)

        with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/profile_dict.json', 'w') as save_as_dict:
            json.dump(profile_dict, save_as_dict)

        with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Soft_SE_dict.json', 'w') as save_as_dict:
            json.dump(Soft_SE_dict, save_as_dict)

        with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/SE_dict.json', 'w') as save_as_dict:
            json.dump(SE_dict, save_as_dict)

    return profile_dict, Soft_SE_dict, SE_dict


def make_prediction(side, FF_dictionary, score_symptom_metric, fixed_symptoms_dict=0):

    ''' Predict symptom-profile improvement for a given activation profile based on the target activation profiles
        fixed symptom dictionary is only needed for weight optimization'''

    profile_dict, Soft_SE_dict, SE_dict = create_NB_dictionaries(side, FF_dictionary, disease='spontaneous human combustion')

    # load fixed weights (they play role only for network blending, not simple prediction)
    if fixed_symptoms_dict == 0:
        fixed_symptom_weights = []  # placeholder
    else:
        with open(fixed_symptoms_dict, 'r') as fp:
            fixed_symptom_weights = json.load(fp)
        fp.close()

    # we can just check the distance for one activation profile across all simulated fibers
    from NB_outline import load_AP_from_OSSDBS
    activation_profile, Pathways = load_AP_from_OSSDBS(side, inters_as_stim=False)

    # get symptom-wise difference between activation and target profiles
    # we do not need non-fixed symptom distances here
    from Optim_strategies import get_symptom_distances
    [__, symptom_diff, symptom_list] = get_symptom_distances(activation_profile, profile_dict, Soft_SE_dict,
                                               fixed_symptom_weights, Pathways, side, score_symptom_metric=score_symptom_metric)
    # symptom_diff is in the symptom space, not pathway! So it might have a different dimensionality

    from Pathways_Stats import get_current_protocol
    current_protocol = get_current_protocol(side)

    from RoutinesForResults import get_activation_prediction
    get_activation_prediction(current_protocol, activation_profile, Pathways, symptom_diff, profile_dict,
                              Soft_SE_dict, side, plot_results=True, score_symptom_metric=score_symptom_metric)

if __name__ == '__main__':

    # called from MATLAB
    # sys.argv[1] - stim folder
    # sys.argv[2] - side
    # sys.argv[3] - Activation Profile Dictionary based on Fiber Filtering
    # sys.argv[4] - Fixed Symptoms Dictionary
    # sys.argv[5] - Score Symptom Metrix

    os.environ['STIMDIR'] = sys.argv[1]
    # make_prediction(int(sys.argv[2]), sys.argv[3], sys.argv[4])
    make_prediction(int(sys.argv[2]),sys.argv[3], sys.argv[5], sys.argv[4])
