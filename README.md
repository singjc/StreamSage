# StreamSage

StreamSage is a web application for performing proteomics database searching using Sage. For more information about Sage, please refer to the github repository [here](https://github.com/lazear/sage).
This app is based on [OpenMS streamlit template project](https://github.com/OpenMS/streamlit-template).

## Installation

Clone the repository and install the dependencies using the following commands:

```bash
git clone https://github.com/singjc/StreamSage.git
cd StreamSage
pip install -r requirements.txt
```

## Usage

Assuming you are in the StreamSage directory

```bash
streamlit run app.py
```

To run the app locally, add the local flag to the command:

**Note**: Running the app locally allows direct access to the users file system. This means you directly use files from your local machine. 

```bash
streamlit run app.py local
```
## Docker

To run the app using Docker, you can pull the available image from Docker Hub:

```bash
docker pull singjust/streamsage:latest
```

Then run the image:

```bash
docker run --rm singjust/streamsage:latest
```

## References

- [Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale](https://doi.org/10.1021/acs.jproteome.3c00486)


