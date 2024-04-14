#!/usr/bin/env bash
# The installation script for makedeb. This is the script that's shown and gets ran from https://makedeb.org.
set -e

# Handy env vars.
MAKEDEB_RELEASE="${MAKEDEB_RELEASE:-}"
makedeb_url='makedeb.org'
color_normal="$(tput sgr0)"
color_bold="$(tput bold)"
color_green="$(tput setaf 77)"
color_orange="$(tput setaf 214)"
color_blue="$(tput setaf 14)"
color_red="$(tput setaf 202)"
color_purple="$(tput setaf 135)"
noninteractive_mode=0
apt_args=()

# Handy functions.
msg() {
    echo "${color_blue}[>]${color_normal} ${1}"
}

error() {
    echo "${color_red}[!]${color_normal} ${1}"
}

question() {
    echo "${color_purple}[?]${color_normal} ${1}"
}

die_cmd() {
    error "${1}"
    exit 1
}

answered_yes() {
    if [[ "${1}" == "" || "${1,,}" == "y" ]]; then
        return 0
    else
        return 1
    fi
}

# Pre-checks.
if [[ "${UID}" == "0" ]]; then
    die_cmd "This script is not allowed to be run under the root user. Please run as a normal user and try again."
fi

# Program start.
echo "-------------------------"
echo "${color_green}[#]${color_normal} ${color_orange}makedeb Installer${color_normal} ${color_green}[#]${color_normal}"
echo "-------------------------"
echo

if ! echo "${-}" | grep -q i; then
    msg "Running in noninteractive mode."
    noninteractive_mode=1
    export DEBIAN_FRONTEND=noninteractive
    apt_args+=('-y')
fi

msg "Ensuring needed packages are installed..."
if ! sudo apt-get update "${apt_args[@]}"; then
    die_cmd "Failed to update APT cache."
fi

missing_dependencies=()
dpkg-query -W 'wget' > /dev/null 2>&1 || missing_dependencies+=('wget')
dpkg-query -W 'gpg' > /dev/null 2>&1 || missing_dependencies+=('gpg')

if ! ( test -z "${missing_dependencies[*]}" || sudo apt-get install "${apt_args[@]}" --mark-auto "${missing_dependencies[@]}" ); then
    die_cmd "Failed to install needed packages."
fi

echo

if (( "${noninteractive_mode}" )) && [[ "${MAKEDEB_RELEASE:+x}" == '' ]]; then
    error "The script was ran in noninteractive mode, but no makedeb package was specified to install."
    error "Please specify a package to install via the 'MAKEDEB_RELEASE' environment variable."
    die_cmd "Available packages are 'makedeb', 'makedeb-beta', and 'makedeb-alpha'."
elif [[ "${MAKEDEB_RELEASE:+x}" == '' ]]; then
    msg "Multiple releases of makedeb are available for installation."
    msg "Currently, you can install one of 'makedeb', 'makedeb-beta', or"
    msg "'makedeb-alpha'."
    
    while true; do
        read -p "$(question "Which release would you like? ")" MAKEDEB_RELEASE

        if ! echo "${MAKEDEB_RELEASE}" | grep -qE '^makedeb$|^makedeb-beta$|^makedeb-alpha$'; then
            error "Invalid response: ${MAKEDEB_RELEASE}"
            continue
        fi

        break
    done
    
    echo
fi

case "${MAKEDEB_RELEASE}" in
    makedeb|makedeb-alpha|makedeb-beta)
        ;;
    *)
        echo
        error "Invalid \$MAKEDEB_RELEASE: '${MAKEDEB_RELEASE}'"
        exit 1 ;;
esac

msg "Setting up makedeb APT repository..."
if ! wget -qO - "https://proget.${makedeb_url}/debian-feeds/makedeb.pub" | gpg --dearmor | sudo tee /usr/share/keyrings/makedeb-archive-keyring.gpg 1> /dev/null; then
    die_cmd "Failed to set up makedeb APT repository."
fi
echo "deb [signed-by=/usr/share/keyrings/makedeb-archive-keyring.gpg arch=all] https://proget.${makedeb_url} makedeb main" | sudo tee /etc/apt/sources.list.d/makedeb.list 1> /dev/null

msg "Updating APT cache..."
if ! sudo apt-get update "${apt_args[@]}"; then
    die_cmd "Failed to update APT cache."
fi

echo
msg "Installing '${MAKEDEB_RELEASE}'..."
if ! sudo apt-get install "${apt_args[@]}" -- "${MAKEDEB_RELEASE}"; then
    die_cmd "Failed to install package."
fi

msg "Finished! If you need help of any kind, feel free to reach out:"
echo
msg "${color_bold}makedeb Homepage:${color_normal}            https://${makedeb_url}"
msg "${color_bold}makedeb Package Repository:${color_normal}  https://mpr.${makedeb_url}"
msg "${color_bold}makedeb Documentation:${color_normal}       https://docs.${makedeb_url}"
msg "${color_bold}makedeb Support:${color_normal}             https://docs.${makedeb_url}/support/obtaining-support"
echo
msg "Enjoy makedeb!"

# vim: set sw=4 expandtab:
