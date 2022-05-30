def check_pvt_charge_temp(htf_temp, temp_lift, top_layer_temp):
    return htf_temp + temp_lift >= top_layer_temp


def charge_storage():
    print("STORAGE CHARGING")


def main():

    templifts = [5, 5, 5, 5, 5, 5]
    bottom_layer_temp = 10
    htf_temp = 10
    top_layer_temp = 25

    for templift in templifts:
        decision_var = check_pvt_charge_temp(htf_temp, templift, top_layer_temp)
        if decision_var:
            charge_storage()
            htf_temp = bottom_layer_temp
            print(htf_temp)
        else:
            htf_temp += templift
            print(htf_temp)


main()
