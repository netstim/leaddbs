import json
import os
import sys

def create_NB_dictionaries(side, FF_dictionary, disease='spontaneous human combustion'):

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


def make_prediction(side, FF_dictionary, fixed_symptoms_dict):

    profile_dict, Soft_SE_dict, SE_dict = create_NB_dictionaries(os.environ['STIMDIR'], side, FF_dictionary, disease='spontaneous human combustion')

    # load fixed weights
    with open(fixed_symptoms_dict, 'r') as fp:
        fixed_symptom_weights = json.load(fp)
    fp.close()

    # we can just check the distance for one activation profile across all simulated fibers
    from NB_outline import load_AP_from_LeadDBS
    activation_profile, Pathways = load_AP_from_LeadDBS(side, inters_as_stim=False)

    # get symptom-wise difference between activation and target profiles
    from Optim_strategies import get_symptom_distances
    [__, symptom_diff] = get_symptom_distances(activation_profile, profile_dict, Soft_SE_dict,
                                               fixed_symptom_weights, Pathways, side, score_symptom_metric='Canberra')
    # symptom_diff is in the symptom space, not pathway! So it might have a smaller dimensionality

    from Pathways_Stats import get_current_protocol
    current_protocol = get_current_protocol(side)

    from RoutinesForResults import get_activation_prediction
    get_activation_prediction(current_protocol, activation_profile, Pathways, symptom_diff, profile_dict,
                              Soft_SE_dict, side, score_symptom_metric='Canberra')

if __name__ == '__main__':

    # called from MATLAB
    # sys.argv[1] - stim folder
    # sys.argv[2] - side
    # sys.argv[3] - Activation Profile Dictionary based on Fiber Filtering
    # sys.argv[4] - Fixed Symptoms

    os.environ['STIMDIR'] = sys.argv[1]
    make_prediction(sys.argv[2], sys.argv[3], sys.argv[4])
